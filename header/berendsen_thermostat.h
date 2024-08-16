//
// Created by asser on 7/6/2024.
//

#ifndef HPC_PROJECT_BERENDSEN_THERMOSTAT_H
#define HPC_PROJECT_BERENDSEN_THERMOSTAT_H

#include "atoms.h"

void berendsen_thermostat(Atoms &atoms, double target_temperature, double timestep, double relaxation_time);



#endif // HPC_PROJECT_BERENDSEN_THERMOSTAT_H
