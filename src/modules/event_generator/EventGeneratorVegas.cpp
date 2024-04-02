/*
 * EventGeneratorVegas.cpp
 *
 *  Created on: Apr 2, 2024
 *      Author: Laurent Forthomme (AGH)
 */

#include <ElementaryUtils/logger/CustomException.h>
#include <ElementaryUtils/parameters/GenericType.h>
#include <ElementaryUtils/string_utils/Formatter.h>
#include <TFile.h>
#include <TRandom3.h>
#include <partons/BaseObjectRegistry.h>
#include <partons/Partons.h>
#include <partons/ServiceObjectRegistry.h>
#include <partons/services/hash_sum/CryptographicHashService.h>
#include <sys/stat.h>

#include "../../../include/beans/other/EventGeneratorInterface.h"
#include "../../../include/modules/event_generator/EventGeneratorVegas.h"

namespace EPIC {

  const std::string EventGeneratorVegas::PARAMETER_NAME_NUM_FCT_CALLS = "numFunctionCalls";
  const std::string EventGeneratorVegas::PARAMETER_NAME_WARMUP_CALLS = "numWarmupCalls";
  const std::string EventGeneratorVegas::PARAMETER_NAME_ITERATIONS = "numIterations";
  const std::string EventGeneratorVegas::PARAMETER_NAME_CHISQ_CUT = "chiSquareCut";
  const std::string EventGeneratorVegas::PARAMETER_NAME_ALPHA = "alpha";
  const std::string EventGeneratorVegas::PARAMETER_NAME_VERBOSITY = "verbosity";
  const std::string EventGeneratorVegas::PARAMETER_NAME_MODE = "mode";

  const unsigned int EventGeneratorVegas::classId =
      PARTONS::BaseObjectRegistry::getInstance()->registerBaseObject(new EventGeneratorVegas("EventGeneratorVegas"));

  EventGeneratorInterface *EventGeneratorVegas::m_pEventGeneratorInterface = nullptr;

  EventGeneratorVegas::EventGeneratorVegas(const std::string &className) : EventGeneratorModule(className) {}

  EventGeneratorVegas::EventGeneratorVegas(const EventGeneratorVegas &other) : EventGeneratorModule(other) {
    if (other.m_vegas_state)
      warn(__func__, "Not able to copy Vegas state object, you need to run initialization for new object");

    m_num_function_calls = other.m_num_function_calls;
    m_num_warmup_calls = other.m_num_warmup_calls;
    m_max_iterations = other.m_max_iterations;
    m_chisq_cut = other.m_chisq_cut;
    m_alpha = other.m_alpha;
    m_verbosity = other.m_verbosity;
    m_mode = other.m_mode;
  }

  EventGeneratorVegas::~EventGeneratorVegas() {}

  EventGeneratorVegas *EventGeneratorVegas::clone() const { return new EventGeneratorVegas(*this); }

  void EventGeneratorVegas::configure(const ElemUtils::Parameters &parameters) {
    EventGeneratorModule::configure(parameters);
    gsl_rng_env_setup();
    const gsl_rng_type *T;
    T = gsl_rng_default;
    m_rnd_gen = gsl_rng_alloc(T);

    if (parameters.isAvailable(EventGeneratorVegas::PARAMETER_NAME_NUM_FCT_CALLS)) {
      m_num_function_calls = parameters.getLastAvailable().toUInt();
      info(__func__,
           ElemUtils::Formatter() << "Parameter " << EventGeneratorVegas::PARAMETER_NAME_NUM_FCT_CALLS << " changed to "
                                  << m_num_function_calls);
    }

    if (parameters.isAvailable(EventGeneratorVegas::PARAMETER_NAME_WARMUP_CALLS)) {
      m_num_warmup_calls = parameters.getLastAvailable().toUInt();
      info(__func__,
           ElemUtils::Formatter() << "Parameter " << EventGeneratorVegas::PARAMETER_NAME_WARMUP_CALLS << " changed to "
                                  << m_num_warmup_calls);
    }

    if (parameters.isAvailable(EventGeneratorVegas::PARAMETER_NAME_ITERATIONS)) {
      m_max_iterations = parameters.getLastAvailable().toUInt();
      info(__func__,
           ElemUtils::Formatter() << "Parameter " << EventGeneratorVegas::PARAMETER_NAME_ITERATIONS << " changed to "
                                  << m_max_iterations);
    }

    if (parameters.isAvailable(EventGeneratorVegas::PARAMETER_NAME_CHISQ_CUT)) {
      m_chisq_cut = parameters.getLastAvailable().toDouble();
      info(__func__,
           ElemUtils::Formatter() << "Parameter " << EventGeneratorVegas::PARAMETER_NAME_CHISQ_CUT << " changed to "
                                  << m_chisq_cut);
    }

    if (parameters.isAvailable(EventGeneratorVegas::PARAMETER_NAME_ALPHA)) {
      m_alpha = parameters.getLastAvailable().toDouble();
      info(__func__,
           ElemUtils::Formatter() << "Parameter " << EventGeneratorVegas::PARAMETER_NAME_ALPHA << " changed to "
                                  << m_alpha);
    }

    if (parameters.isAvailable(EventGeneratorVegas::PARAMETER_NAME_VERBOSITY)) {
      m_verbosity = parameters.getLastAvailable().toUInt();
      info(__func__,
           ElemUtils::Formatter() << "Parameter " << EventGeneratorVegas::PARAMETER_NAME_VERBOSITY << " changed to "
                                  << m_verbosity);
    }

    if (parameters.isAvailable(EventGeneratorVegas::PARAMETER_NAME_MODE)) {
      m_mode = parameters.getLastAvailable().toUInt();
      info(__func__,
           ElemUtils::Formatter() << "Parameter " << EventGeneratorVegas::PARAMETER_NAME_MODE << " changed to "
                                  << m_mode);
    }
  }

  void EventGeneratorVegas::initialise(const std::vector<KinematicRange> &kinematicRanges,
                                       const EventGeneratorInterface &service) {
    //scenario
    size_t scenario = 0;

    if (!m_initStatePath.empty()) {
      struct stat buffer;
      if (!stat(m_initStatePath.c_str(), &buffer) == 0)  // file does not exist
        scenario = 1;
      else  // file exists
        scenario = 2;
    }

    m_kinematicRanges = kinematicRanges;
    m_pEventGeneratorInterface = const_cast<EventGeneratorInterface *>(&service);

    // initialize
    if (scenario == 0 || scenario == 1) {
      info(__func__, "Creating new Vegas integrator object");

      // integrator parameters
      m_vegas_state.reset(gsl_monte_vegas_alloc(m_kinematicRanges.size()));
      gsl_monte_vegas_params_get(m_vegas_state.get(), &m_vegas_params);
      m_vegas_params.iterations = m_max_iterations;
      m_vegas_params.alpha = m_alpha;
      m_vegas_params.verbose = m_verbosity;
      m_vegas_params.mode = m_mode;
      gsl_monte_vegas_params_set(m_vegas_state.get(), &m_vegas_params);

      if (m_pEventGeneratorInterface == nullptr)
        throw ElemUtils::CustomException(getClassName(), __func__, "Pointer to EventGenerator is null");
      if (m_kinematicRanges.size() != kinematicRanges.size())
        throw ElemUtils::CustomException(
            getClassName(), __func__, "Size of vector containing  kinematic ranges different than nDim");

      // function definition
      m_function.f = integrand;
      m_function.dim = m_kinematicRanges.size();
      m_function.params = nullptr;

      // phase space definition
      m_lowx.resize(m_kinematicRanges.size());
      m_highx.resize(m_kinematicRanges.size());
      for (size_t i = 0; i < m_kinematicRanges.size(); ++i) {
        m_lowx[i] = m_kinematicRanges.at(i).getMin();
        m_highx[i] = m_kinematicRanges.at(i).getMax();
      }

      // compute the integral
      warmup();
      size_t it_chisq = 0;
      double result = 0., abserr = 0.;
      do {  // main integration run
        iterate(m_num_function_calls / m_max_iterations, result, abserr);
        ++it_chisq;
        info(__func__,
             ElemUtils::Formatter() << "Iteration #" << it_chisq << ": average=" << result << ", sigma=" << abserr
                                    << ", chi2=" << gsl_monte_vegas_chisq(m_vegas_state.get()) << ".");
      } while (std::fabs(gsl_monte_vegas_chisq(m_vegas_state.get()) - 1.) > m_chisq_cut - 1.);
      m_integral = std::make_pair(result, abserr);
    }
  }

  void EventGeneratorVegas::warmup() {
    double result = 0., abserr = 0.;
    iterate(m_num_warmup_calls, result, abserr);
    info(__func__,
         ElemUtils::Formatter() << "Vegas finished warmup with " << m_num_warmup_calls
                                << " calls. Initial value of the integral: " << result << " +/- " << abserr << ".");
    m_integral = std::make_pair(result, abserr);
  }

  void EventGeneratorVegas::iterate(size_t num_calls, double &result, double &abserr) {
    if (int res = gsl_monte_vegas_integrate(&m_function,
                                            &m_lowx[0],
                                            &m_highx[0],
                                            m_kinematicRanges.size(),
                                            num_calls,
                                            m_rnd_gen,
                                            m_vegas_state.get(),
                                            &result,
                                            &abserr);
        res != GSL_SUCCESS)
      throw ElemUtils::CustomException(getClassName(),
                                       __func__,
                                       ElemUtils::Formatter() << "Failed to iterate with " << num_calls
                                                              << " calls. GSL error: " << gsl_strerror(res) << ".");
  }

  double EventGeneratorVegas::integrand(double *x, size_t dim, void *params) {
    std::vector<double> v(x, x + dim);
    return m_pEventGeneratorInterface->getEventDistribution(v);
  }

  std::pair<std::vector<double>, double> EventGeneratorVegas::generateEvent() {
    throw ElemUtils::CustomException(getClassName(), __func__, "Module is not yet ready for event generation");
  }

  std::pair<double, double> EventGeneratorVegas::getIntegral() { return m_integral; }

} /* namespace EPIC */
