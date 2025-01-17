/*
 * WriterHepMC3.cpp
 *
 *  Created on: Feb 11, 2021
 *      Author: Pawel Sznajder (NCBJ) 
 */

#include "../../../include/modules/writer/WriterHepMC3.h"

#include <ElementaryUtils/logger/CustomException.h>
#include <ElementaryUtils/parameters/GenericType.h>
#include <ElementaryUtils/string_utils/Formatter.h>
#include <HepMC3/Attribute.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/GenCrossSection.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Units.h>
#include <HepMC3/Writer.h>
#include <HepMC3/WriterAscii.h>
#include <HepMC3/WriterHEPEVT.h>
#if WITH_HEPMC3ROOTIO
	#include <HepMC3/WriterRoot.h>
#endif
#include <partons/BaseObjectRegistry.h>
#include <stddef.h>
#include <TLorentzVector.h>
#include <iterator>
#include <map>
#include <utility>

#include "../../../include/beans/containers/GenerationInformation.h"
#include "../../../include/beans/physics/Particle.h"
#include "../../../include/beans/physics/Vertex.h"
#include "../../../include/beans/types/EventAttributeType.h"
#include "../../../include/beans/types/ParticleCodeType.h"

namespace EPIC {

const std::string WriterHepMC3::PARAMETER_NAME_HEPMC3_WRITER_TYPE =
		"HepMC3_writer_type";

const unsigned int WriterHepMC3::classId =
		PARTONS::BaseObjectRegistry::getInstance()->registerBaseObject(
				new WriterHepMC3("WriterHepMC3"));

WriterHepMC3::WriterHepMC3(const std::string &className) :
		WriterModule(className), m_writerHepMC3Type(
				WriterHepMC3Type::UNDEFINED), m_writerHepMC3(nullptr), m_runInfoHepMC3(nullptr) {
}

WriterHepMC3::WriterHepMC3(const WriterHepMC3 &other) :
		WriterModule(other), m_writerHepMC3Type(other.m_writerHepMC3Type), m_writerHepMC3(
				other.m_writerHepMC3), m_runInfoHepMC3(other.m_runInfoHepMC3) {
	;
}

WriterHepMC3::~WriterHepMC3() {
}

WriterHepMC3 *WriterHepMC3::clone() const {
	return new WriterHepMC3(*this);
}

void WriterHepMC3::configure(const ElemUtils::Parameters &parameters) {

	// run for mother
	WriterModule::configure(parameters);

	// look for name
	if (parameters.isAvailable(
			WriterHepMC3::PARAMETER_NAME_HEPMC3_WRITER_TYPE)) {

		// get
		m_writerHepMC3Type = WriterHepMC3Type::fromString(
				parameters.getLastAvailable().getString());

		// print status
		info(__func__,
				ElemUtils::Formatter() << "HepMC3 writer type set to: "
						<< WriterHepMC3Type(m_writerHepMC3Type).toString());
	}
}

void WriterHepMC3::open() {

	// check path
	if (m_path.empty()) {
		throw ElemUtils::CustomException(getClassName(), __func__,
				"Path not set");
	}

	// close if was open
	if (m_writerHepMC3 != nullptr) {
		close();
	}

	// switch
	switch (m_writerHepMC3Type) {

	case WriterHepMC3Type::ASCII: {

		m_writerHepMC3 = std::make_shared<HepMC3::WriterAscii>(m_path);
		break;
	}

#if WITH_HEPMC3ROOTIO
	case WriterHepMC3Type::ROOT: {

		m_writerHepMC3 = std::make_shared<HepMC3::WriterRoot>(m_path);
		break;
	}
#endif

	case WriterHepMC3Type::HEPEVT: {

		m_writerHepMC3 = std::make_shared<HepMC3::WriterHEPEVT>(m_path);
		break;
	}

	default: {
		throw ElemUtils::CustomException(getClassName(), __func__,
				ElemUtils::Formatter() << "HepMC3 writer type is: "
						<< WriterHepMC3Type(m_writerHepMC3Type).toString());
	}
	}

	//delete if already set
	if (m_runInfoHepMC3 != nullptr) {
		m_runInfoHepMC3.reset(new HepMC3::GenRunInfo());
	}else{
		m_runInfoHepMC3 = std::make_shared<HepMC3::GenRunInfo>();
	}

	//set weight
	std::vector<std::string> weights(1, "weight_default");
	m_runInfoHepMC3->set_weight_names(weights);

	//assign
	m_writerHepMC3->set_run_info(m_runInfoHepMC3);
}

void WriterHepMC3::saveGenerationInformation(
		const GenerationInformation& generationInformation) {

	if (m_writerHepMC3 == nullptr) {

		warn(__func__, "HepMC3 writer is not open");
		return;
	}

	HepMC3::GenRunInfo::ToolInfo toolInfo;

	toolInfo.name = generationInformation.getGeneratorName();
	toolInfo.version = generationInformation.getGeneratorVersion();
	toolInfo.description = generationInformation.getDescription();

	m_runInfoHepMC3->tools().push_back(toolInfo);

	std::shared_ptr<HepMC3::StringAttribute> attributeServiceName =
			std::make_shared<HepMC3::StringAttribute>(
					generationInformation.getServiceName());
	m_runInfoHepMC3->add_attribute("service_name", attributeServiceName);

	std::shared_ptr<HepMC3::UIntAttribute> attributeNEvents = std::make_shared<
			HepMC3::UIntAttribute>(generationInformation.getNEvents());
	m_runInfoHepMC3->add_attribute("generated_events_number", attributeNEvents);

	std::shared_ptr<HepMC3::DoubleAttribute> attributeICS_value =
			std::make_shared<HepMC3::DoubleAttribute>(
					generationInformation.getIntegratedCrossSection().first);
	m_runInfoHepMC3->add_attribute("integrated_cross_section_value",
			attributeICS_value);

	std::shared_ptr<HepMC3::DoubleAttribute> attributeICS_uncertainty =
			std::make_shared<HepMC3::DoubleAttribute>(
					generationInformation.getIntegratedCrossSection().second);
	m_runInfoHepMC3->add_attribute("integrated_cross_section_uncertainty",
			attributeICS_uncertainty);

	std::shared_ptr<HepMC3::StringAttribute> attributeGenerationDate =
			std::make_shared<HepMC3::StringAttribute>(
					generationInformation.getGenerationDate());
	m_runInfoHepMC3->add_attribute("generation_date", attributeGenerationDate);

	for (std::vector<std::pair<std::string, std::string> >::const_iterator it =
			generationInformation.getAdditionalInfo().begin();
			it != generationInformation.getAdditionalInfo().end(); it++) {

		std::shared_ptr<HepMC3::StringAttribute> attributeAdditional =
				std::make_shared<HepMC3::StringAttribute>(it->second);
		m_runInfoHepMC3->add_attribute(it->first, attributeAdditional);
	}

	switch (m_writerHepMC3Type) {

	case WriterHepMC3Type::ASCII: {

		std::static_pointer_cast<HepMC3::WriterAscii>(m_writerHepMC3)->write_run_info();
		break;
	}

#if WITH_HEPMC3ROOTIO
	case WriterHepMC3Type::ROOT: {

		std::static_pointer_cast<HepMC3::WriterRoot>(m_writerHepMC3)->write_run_info();
		break;
	}
#endif

	case WriterHepMC3Type::HEPEVT: {

		warn(__func__,
				ElemUtils::Formatter()
						<< "Saving generation information not available for writer type "
						<< WriterHepMC3Type(m_writerHepMC3Type).toString());
		break;
	}

	default: {
		throw ElemUtils::CustomException(getClassName(), __func__,
				ElemUtils::Formatter() << "HepMC3 writer type is: "
						<< WriterHepMC3Type(m_writerHepMC3Type).toString());
	}
	}
}

void WriterHepMC3::close() {

	if (m_writerHepMC3 == nullptr) {

		warn(__func__, "HepMC3 writer is not open");
		return;
	}

	m_writerHepMC3->close();
}

void WriterHepMC3::write(const Event &event) {

	// check if open
	if (m_writerHepMC3 == nullptr) {

		throw ElemUtils::CustomException(getClassName(), __func__,
				"HepMC3 writer is not open");
	}

	// references
	const std::vector<std::pair<ParticleCodeType::Type, std::shared_ptr<Particle> > > &eventParticles =
			event.getParticles();
	const std::vector<std::shared_ptr<Vertex> > &eventVertices = event.getVertices();
	const std::map<EventAttributeType::Type, double> &eventAttributes =
			event.getAttributes();

	// event
	HepMC3::GenEvent eventHepMC3;

	// set units
	eventHepMC3.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);

	//set run info
	eventHepMC3.set_run_info(m_runInfoHepMC3);

	// attributes
	for (std::map<EventAttributeType::Type, double>::const_iterator it =
			eventAttributes.begin(); it != eventAttributes.end(); it++) {

		// cross section
		if (it->first == EventAttributeType::WEIGHT) {
			std::shared_ptr<HepMC3::GenCrossSection> cross_section =
					std::make_shared<HepMC3::GenCrossSection>();
			cross_section->set_cross_section(it->second, 0.);
			eventHepMC3.add_attribute("GenCrossSection", cross_section);
		}

		// set number
		if (it->first == EventAttributeType::ID) {
			eventHepMC3.set_event_number(int(it->second));
		}
	}

	// pairs
	std::vector<
			std::pair<std::shared_ptr<Particle>,
					std::shared_ptr<HepMC3::GenParticle>>> mc3pairs(
			eventParticles.size());

	// particles
	for (std::vector<std::pair<ParticleCodeType::Type, std::shared_ptr<Particle> > >::const_iterator it =
			eventParticles.begin(); it != eventParticles.end(); it++) {

		// create pair
		mc3pairs.at(size_t(it - eventParticles.begin())) = std::make_pair(it->second,
				std::make_shared<HepMC3::GenParticle>(
						HepMC3::FourVector(it->second->getFourMomentum().Px(),
								it->second->getFourMomentum().Py(),
								it->second->getFourMomentum().Pz(),
								it->second->getFourMomentum().E()),
						(int) (it->second->getType()),
						getParticleCode(it->first)));
	}

	// vertices
	for (std::vector<std::shared_ptr<Vertex> >::const_iterator it0 = eventVertices.begin();
			it0 != eventVertices.end(); it0++) {

		// create vertex
		std::shared_ptr<HepMC3::GenVertex> mc3Vertex = std::make_shared<
				HepMC3::GenVertex>();

		// in
		const std::vector<std::shared_ptr<Particle>> &eventVertexParticlesIn =
				(*it0)->getParticlesIn();

		for (std::vector<std::shared_ptr<Particle>>::const_iterator it1 =
				eventVertexParticlesIn.begin();
				it1 != eventVertexParticlesIn.end(); it1++) {
			for (std::vector<
					std::pair<std::shared_ptr<Particle>,
							std::shared_ptr<HepMC3::GenParticle>>>::const_iterator it2 =
					mc3pairs.begin(); it2 != mc3pairs.end(); it2++) {
				if ((*it1) == it2->first) {

					mc3Vertex->add_particle_in(it2->second);
					break;
				}
			}
		}

		// out
		const std::vector<std::shared_ptr<Particle>> &eventVertexParticlesOut =
				(*it0)->getParticlesOut();

		for (std::vector<std::shared_ptr<Particle>>::const_iterator it1 =
				eventVertexParticlesOut.begin();
				it1 != eventVertexParticlesOut.end(); it1++) {
			for (std::vector<
					std::pair<std::shared_ptr<Particle>,
							std::shared_ptr<HepMC3::GenParticle>>>::const_iterator it2 =
					mc3pairs.begin(); it2 != mc3pairs.end(); it2++) {
				if ((*it1) == it2->first) {

					mc3Vertex->add_particle_out(it2->second);
					break;
				}
			}
		}

		// add
		eventHepMC3.add_vertex(mc3Vertex);
	}

	// write
	m_writerHepMC3->write_event(eventHepMC3);
}

void WriterHepMC3::write(const std::vector<Event> &events) {

	for (std::vector<Event>::const_iterator it = events.begin();
			it != events.end(); it++) {
		write(*it);
	}
}

int WriterHepMC3::getParticleCode(ParticleCodeType::Type type) const {

    switch (type) {

    case ParticleCodeType::UNDECAYED:
        return 1;
        break;
    case ParticleCodeType::DECAYED:
        return 2;
        break;
    case ParticleCodeType::DOCUMENTATION:
        return 3;
        break;
    case ParticleCodeType::BEAM:
        return 4;
        break;
    case ParticleCodeType::SCATTERED:
        return 1;
        break;
    case ParticleCodeType::VIRTUAL:
        return 13;
        break;

    default:
        throw ElemUtils::CustomException(getClassName(), __func__,
                ElemUtils::Formatter() << "Conversion undefined for type: "
                        << ParticleCodeType(type).toString());
    }
}

WriterHepMC3Type::Type WriterHepMC3::getWriterHepMC3Type() const {
	return m_writerHepMC3Type;
}

void WriterHepMC3::setWriterHepMC3Type(
		WriterHepMC3Type::Type writerHepMc3Type) {

	if (m_writerHepMC3 != nullptr) {

		warn(__func__,
				ElemUtils::Formatter()
						<< "HepMC3 writer already is open with type:"
						<< WriterHepMC3Type(m_writerHepMC3Type).toString()
						<< ", change to "
						<< WriterHepMC3Type(writerHepMc3Type).toString()
						<< " ignored");
		return;
	}

	m_writerHepMC3Type = writerHepMc3Type;
}

} /* namespace EPIC */
