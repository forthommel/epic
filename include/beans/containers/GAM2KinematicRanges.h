/*
 * GAM2KinematicRanges.h
 *
 *  Created on: Feb 9, 2021
 *      Author: Pawel Sznajder (NCBJ)
 */

#ifndef BEANS_GAM2KINEMATICRANGES_H_
#define BEANS_GAM2KINEMATICRANGES_H_

#include <string>

#include "KinematicRange.h"

namespace EPIC {
class GAM2Kinematic;
class ExperimentalConditions;
} /* namespace EPIC */

namespace EPIC {

/**
 * @class GAM2KinematicRanges
 *
 * @brief Container to store single kinematic ranges for GAM2.
 *
 * This class acts as a container to store kinematic range for GAM2.
 */
class GAM2KinematicRanges: public PARTONS::BaseObject {

public:

    static const std::string GAM2_KINEMATIC_RANGE_Y; ///< Key to set y kinematic range.
    static const std::string GAM2_KINEMATIC_RANGE_Q2; ///< Key to set Q2 kinematic range.
    static const std::string GAM2_KINEMATIC_RANGE_T; ///< Key to set t kinematic range.
    static const std::string GAM2_KINEMATIC_RANGE_UPRIM; ///< Key to set u' kinematic range.
    static const std::string GAM2_KINEMATIC_RANGE_MGG2; ///< Key to set Mgg2 kinematic range.
    static const std::string GAM2_KINEMATIC_RANGE_PHI; ///< Key to set phi kinematic range.

    /**
     * Default constructor.
     */
    GAM2KinematicRanges();

    /**
     * Assignment constructor.
     */
    GAM2KinematicRanges(const KinematicRange &rangeY,
            const KinematicRange &rangeQ2, const KinematicRange &rangeT,
            const KinematicRange &rangeUPrim, const KinematicRange &rangeMgg2,
            const KinematicRange &rangePhi);

    /**
     * Copy constructor.
     */
    GAM2KinematicRanges(const GAM2KinematicRanges &other);

    /**
     * Destructor.
     */
    virtual ~GAM2KinematicRanges();

    virtual std::string toString() const;

    bool inRange(const ExperimentalConditions& experimentalConditions,
            const GAM2Kinematic& obsKin) const;

    /**
     * Build from task.
     */
    static GAM2KinematicRanges getGAM2KinematicRangesfromTask(const MonteCarloTask &task);

    //********************************************************
    //*** SETTERS AND GETTERS ********************************
    //********************************************************

    /**
     * Get kinematic range y.
     */
    const KinematicRange& getRangeY() const;

    /**
     * Set kinematic range y.
     */
    void setRangeY(const KinematicRange& rangeY);

    /**
     * Get kinematic range Q2.
     */
    const KinematicRange& getRangeQ2() const;

    /**
     * Set kinematic range Q2.
     */
    void setRangeQ2(const KinematicRange& rangeQ2);

    /**
     * Get kinematic range t.
     */
    const KinematicRange& getRangeT() const;

    /**
     * Set kinematic range t.
     */
    void setRangeT(const KinematicRange& rangeT);

    /**
     * Get kinematic range u'.
     */
    const KinematicRange& getRangeUPrim() const;

    /**
     * Set kinematic range u'.
     */
    void setRangeUPrim(const KinematicRange& rangeUPrim);

    /**
     * Get kinematic range Mgg2.
     */
    const KinematicRange& getRangeMgg2() const;

    /**
     * Set kinematic range Mgg2.
     */
    void setRangeMgg2(const KinematicRange& rangeMgg2);

    /**
     * Get kinematic range phi.
     */
    const KinematicRange& getRangePhi() const;

    /**
     * Set kinematic range phi.
     */
    void setRangePhi(const KinematicRange& rangePhi);

private:

    KinematicRange m_rangeY;   ///< Kinematic range y.
    KinematicRange m_rangeQ2;   ///< Kinematic range Q2.
    KinematicRange m_rangeT;   ///< Kinematic range t.
    KinematicRange m_rangeUPrim;  ///< Kinematic range u'.
    KinematicRange m_rangeMgg2; ///< Kinematic range Mgg2.
    KinematicRange m_rangePhi;   ///< Kinematic range phi.
};

} /* namespace EPIC */

#endif /* BEANS_GAM2KINEMATICRANGES_H_ */
