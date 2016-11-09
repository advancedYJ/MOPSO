#include "MOPSO.h"

void initializeParticle(Particle *particle){
    for (int j=0; j<nPop; j++){
        //  initialize Velocity
        int lenV =particle[j].sizeOfPosition;
        particle[j].Velocity = new double[lenV];
        for (int i=0; i<lenV; i++)
            particle[j].Velocity[i] = 0;
        particle[j].sizeOfVelocity =lenV;

        //  old_position
        particle[j].old_position = new double[lenV];

        //  Best
        particle[j].Best.Position = new double [lenV];
        cpyDoubleArray(particle[j].Best.Position,particle[j].Position,lenV);
        particle[j].Best.Cost = new double[numObjective];
        cpyDoubleArray(particle[j].Best.Cost, particle[j].Cost, numObjective);
    }
}
