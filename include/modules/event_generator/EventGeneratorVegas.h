/*
 * EventGeneratorVegas.h
 *
 *  Created on: Apr 2, 2024
 *      Author: Laurent Forthomme (AGH)
 */

#ifndef MODULES_EVENT_GENERATOR_EVENTGENERATORVEGAS_H_
#define MODULES_EVENT_GENERATOR_EVENTGENERATORVEGAS_H_

#include <ElementaryUtils/parameters/Parameters.h>
#include <Rtypes.h>
#include <gsl/gsl_monte_vegas.h>
#include <stddef.h>

#include <string>
#include <utility>
#include <vector>

#include "../../beans/containers/KinematicRange.h"
#include "EventGeneratorModule.h"

namespace EPIC {
  /**
   * @class EventGeneratorVegas
   *
   * @brief Event generator based on Vegas IntegratorMultiDim.
   */
  class EventGeneratorVegas : public EventGeneratorModule {
  public:
    static const std::string PARAMETER_NAME_NUM_FCT_CALLS;  ///< Key to set EventGeneratorVegas::m_num_function_calls.
    static const std::string PARAMETER_NAME_WARMUP_CALLS;   ///< Key to set EventGeneratorVegas::m_num_warmup_calls.
    static const std::string PARAMETER_NAME_ITERATIONS;     ///< Key to set EventGeneratorVegas::m_iterations.
    static const std::string PARAMETER_NAME_CHISQ_CUT;      ///< Key to set EventGeneratorVegas::m_chisq_cut.
    static const std::string PARAMETER_NAME_ALPHA;          ///< Key to set EventGeneratorVegas::m_alpha.
    static const std::string PARAMETER_NAME_VERBOSITY;      ///< Key to set EventGeneratorVegas::m_verbosity.
    static const std::string PARAMETER_NAME_MODE;           ///< Key to set EventGeneratorVegas::m_mode.

    static const unsigned int classId;  ///< Unique ID to automatically register the class in the registry.

    EventGeneratorVegas(const std::string &className);      ///< Default constructor
    EventGeneratorVegas(const EventGeneratorVegas &other);  ///< Copy constructor
    virtual ~EventGeneratorVegas();                         ///< Destructor

    virtual EventGeneratorVegas *clone() const;
    virtual void configure(const ElemUtils::Parameters &parameters);

    virtual void initialise(const std::vector<KinematicRange> &kinematicRanges, const EventGeneratorInterface &service);
    virtual std::pair<std::vector<double>, double> generateEvent();
    virtual std::pair<double, double> getIntegral();

  private:
    static double integrand(double *, size_t, void *);
    void warmup();

    gsl_rng *m_rnd_gen{nullptr};
    gsl_monte_vegas_params m_vegas_params;
    /// Trivial deleter for the Vegas integrator
    struct gsl_monte_vegas_deleter {
      inline void operator()(gsl_monte_vegas_state *state) { gsl_monte_vegas_free(state); }
    };
    /// Vegas integrator state for integration (optional) and/or "treated" event generation
    std::unique_ptr<gsl_monte_vegas_state, gsl_monte_vegas_deleter> m_vegas_state;
    gsl_monte_function m_function;

    /// Pointer used during initialization of Vegas.
    static EventGeneratorInterface *m_pEventGeneratorInterface;
    std::vector<KinematicRange> m_kinematicRanges;  ///< Kinematic ranges set during initialization.

    size_t m_num_function_calls{100'000};
    size_t m_num_warmup_calls{25'000};
    size_t m_iterations{10};  ///< Number of iterations in standard run
    double m_chisq_cut{1.5};
    double m_alpha{1.25};
    int m_verbosity{-1};
    int m_mode{-1};

    std::vector<double> m_lowx, m_highx;

    std::pair<double, double> m_integral{0., 0.};  ///< Integrated cross section + uncertainty
  };
}  // namespace EPIC

#endif
