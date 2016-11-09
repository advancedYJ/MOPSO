//
// Created by ws on 16-7-15.
//

#include "MOPSO.h"

char * inputSeq(){
    // read in the sequence
    char *seq;
    char buffer[256];
    char *filename;
    ifstream inFile;
    filename = catStrStr(inputAddress,"seq.txt");
    openFile(filename, inFile);
    int cntLines = 0;
    while (! inFile.eof()) { inFile.getline(buffer, 256); cntLines++; }
    cntLines--;
    inFile.close();

    openFile(filename, inFile);
    char tmp[cntLines][256];
    int len=0;
    inFile.getline(buffer, 256);                    //  the first line is not needed
    for (int i=0; i<cntLines-1; i++) {
        inFile.getline(tmp[i], 256);
        len += strlen(tmp[i]);
    }

    seq = new char [len];
    strcpy(seq, "");
    for (int i=0; i<cntLines-1; i++)  strcat(seq, tmp[i]);
    inFile.close();
    return seq;

}

// input the coordinaries and angles
void inputParticle(Particle &particle, int filenum, char *seq){
    ifstream inFile;
    char *fileAddress;
    char *filename;
    fileAddress = catStrStr(inputAddress, "particleO_");
    filename = catStrIntStr(fileAddress, filenum+1, ".txt");
    openFile(filename, inFile);
    char buffer[256];
    // get the number of lines of the input file
    int cntLines = 0;
    while (! inFile.eof()) { inFile.getline(buffer, 256); cntLines++; }
    cntLines--;
    particle.numAA = cntLines/4;                                            // the number of AA
    particle.numAtom = particle.numAA *3;                          // the number of atoms ( not add O)
    inFile.close();                                                                         //get how many lines the input has

    // read the file to get the coordinaries
    openFile(filename, inFile);
    particle.addO_origin = new double[3*cntLines];                  //  the x,y,z of every atoms in the particle
    particle.origin = new double[particle.numAA*3 *3];          //  x,y,z, delete the number of the atom O
    int k=0;
    //  read the file by lines
    for (int j=0; j<12 * particle.numAA; j++){
        inFile >> particle.addO_origin[j];
        if (j % 12<9)
            particle.origin[k++] = particle.addO_origin[j];
    }
    particle.sizeOfAddO_origin = particle.numAA*4*3;
    particle.sizeOfOrigin = particle.numAA*3*3;                 //  the totally number of coordinaries
    inFile.close();
    delete filename;

    //  read the init angles
    fileAddress = catStrStr(inputAddress, "phi");
    filename = catStrIntStr(fileAddress, filenum+1, ".txt");
    openFile(filename, inFile);
    particle.Position = new double[3 * particle.numAA-3];     //   why so many 180  ?

    for (int j=0; j< 3*particle.numAA-5; j+=3){          // last : j=3*numAA-6; j+2=3*numAA-4
        inFile >> particle.Position[j] ;
        particle.Position[j+1] = 180;
        inFile >> particle.Position[j+2] ;
    }

    inFile.close();
    particle.sizeOfPosition = 3*particle.numAA-3;

    //  init particle.cost[]
    particle.Cost = new double [numObjective];

    particle.index = filenum;
    particle.seq = seq;

}

void openFile(const char fileDir[10000], ifstream &infile){
    infile.open(fileDir);
    if (!infile) {
        cout << "can't open:::" << fileDir << endl;
        exit(-1);
    }
}

void removeFile(const char * str){
    int k = remove(str);
    if (k==-1) {
        cout << "Fail to delete:::" << str << endl;
        exit(-1);
    }
}

void printPdb(Particle * particle){
    FILE * output;
    //cout << "index= " << particle->index << endl;
    //const char *fileName = catStrIntStr("/home/ws/zzZyj/data/energyFile/tempFile/temp", particle->index, ".pdb");
    const char *fileName = catStrIntStr(energyFileAddress,"tempFile/temp", particle->index, ".pdb");
    output = fopen(fileName, "w");
    if (output == NULL) {
        cout << "Can't write temp.pdb" << endl;
        exit(-1);
    }
    int lenOfAA = particle->numAA;
    for (int i=0; i<lenOfAA; i++){
        char *AA = getRelativeAA(particle->seq[i]);
        int k1 = i*4*3, k2 = i*4*3+3, k3 = i*4*3+6, k4 = i*4*3+9;
        fprintf(output, "ATOM   %4d  N   %3s   %3d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                4*i+1, AA, i+1, particle->addO_origin[k1], particle->addO_origin[k1+1], particle->addO_origin[k1+2]);
        fprintf(output, "ATOM   %4d  CA  %3s   %3d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                4*i+2, AA, i+1, particle->addO_origin[k2], particle->addO_origin[k2+1], particle->addO_origin[k2+2]);
        fprintf(output, "ATOM   %4d  C   %3s   %3d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                4*i+3, AA, i+1, particle->addO_origin[k3], particle->addO_origin[k3+1], particle->addO_origin[k3+2]);
        fprintf(output, "ATOM   %4d  O   %3s   %3d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                4*i+4, AA, i+1, particle->addO_origin[k4], particle->addO_origin[k4+1], particle->addO_origin[k4+2]);
        delete AA;
    }
    fclose(output);
    //system("/home/ws/GL/refine1/scrwl4/Scwrl4 -i /home/ws/zzZyj/data/temp.pdb -o /home/ws/zzZyj/data/temp1.pdb");
    //const char *newFile = catStrIntStr("/home/ws/zzZyj/data/energyFile/tempFile/temp_", particle->index, ".pdb");
    const char *newFile = catStrIntStr(energyFileAddress,"tempFile/temp_", particle->index, ".pdb");
    const char *command = catStrStr(refine_1Address, "scrwl4/Scwrl4 -i ", fileName,
                                    " -o ", newFile);
    system(command);
}

void printParticleCost(Particle * particle, int iterator){
    char *costFile = catStrStr(answerAddress, "cost.txt");
    FILE *outFile;
    outFile = fopen(costFile, "a");
    if (outFile == NULL) {
        cout << "Can't print rep answer" << endl;
        exit(-1);
    }
    fprintf(outFile, "iterator=%d\n", iterator);
    for (int i=0; i<nPop; i++){
        /*fprintf(outFile, "i=%d cost0 =%10.3f cost1 =%10.3f\n",
                i, particle[i].Cost[0], particle[i].Cost[1]);*/
        fprintf(outFile, "i=%d cost0 =%10.3f\n",
                i, particle[i].Cost[0]);
    }
    fclose(outFile);
}

//  fprintf
void printAnswer(myRep & rep){
    myRep ::iterator it1;
    int k=0;
    for (it1=rep.begin(); it1!= rep.end(); it1++){

        char *filename;
        filename = catStrIntStr(answerAddress, "rep", k, ".txt");
        FILE *outFile;
        outFile = fopen(filename, "w");
        if (outFile == NULL) {
            cout << "Can't print rep answer" << endl;
            exit(-1);
        }
        fprintf(outFile, "iterator=%d\n",it1->iterator);
        //fprintf(outFile, "cost0 =%10.3f cost1 =%10.3f\n", it1->Cost[0], it1->Cost[1]);
        fprintf(outFile, "cost0 =%10.3f\n", it1->Cost[0]);
        int len = it1->sizeOfAddO_origin;
        int t=0;
        for (int j=0; j<len; j+=3){
            fprintf(outFile, "%4d ", t);
            switch (t % 4){
                case 0: fprintf(outFile, " N  "); break;
                case 1: fprintf(outFile, " CA "); break;
                case 2: fprintf(outFile, " C  "); break;
                case 3: fprintf(outFile, " O  "); break;
            }
            fprintf(outFile, "%8.3f %8.3f %8.3f\n",
                    it1->addO_origin[j], it1->addO_origin[j+1], it1->addO_origin[j+2]);
            t++;
        }
        k++;
        fclose(outFile);
    }
}

void printTime(int choice, time_t start, time_t finish){
    const char *timeFile = catStrStr(answerAddress, "time.txt");
    ofstream outfile(timeFile, ios::app);
    long duration;
    duration = finish-start;
    switch (choice){
        case 0: outfile << "TOTAL Time: \n"; break;
        case 1: outfile << "Cost 1 : "; break;
        case 2: outfile << "Cost 2 : "; break;
        case 3: outfile << "Transfer Time: "; break;
        case 4: outfile << "Input Time: "; break;
        case 5: outfile << "Adaption Time: "; break;
        case 6: outfile << "IntoRep Time: "; break;
        default: break;
    }
    int seconds = duration % 60;
    duration /= 60;
    int minutes = duration % 60;
    duration /= 60;
    int hours = duration;

    outfile << hours << ":" << minutes << ":" << seconds << endl;
    outfile.close();
}

void createNewFold(){
    const char *mkdirCommand = catStrStr("mkdir ",answerAddress);
    int t= system(mkdirCommand);
    if (t != 0){
        cout << "can't create a new fold for answers" << endl;
        exit(-1);
    }

}

void printParticleCostForDebug(Particle * particle, int iterator){
    char *costFile = catStrStr(draftAddress, "cost.txt");
    FILE *outFile;
    outFile = fopen(costFile, "a");

    if (outFile == NULL) {
        cout << "Can't print rep answer" << endl;
        exit(-1);
    }
    fprintf(outFile, "iterator=%d\n", iterator);
    for (int i=0; i<nPop; i++){
        fprintf(outFile, "i=%d cost0 =%10.3f cost1 =%10.3f\n",
                i, particle[i].Cost[0], particle[i].Cost[1]);
    }
    fclose(outFile);
}

void printAnswerForDebug(myRep & rep){
    myRep ::iterator it1;
    int k=0;
    for (it1=rep.begin(); it1!= rep.end(); it1++){
        char *filename;
        filename = catStrIntStr(draftAddress, "rep", k, ".txt");
        FILE *outFile;
        outFile = fopen(filename, "w");
        if (outFile == NULL) {
            cout << "Can't print rep answer" << endl;
            exit(-1);
        }
        fprintf(outFile, "iterator=%d\n",it1->iterator);
        fprintf(outFile, "cost0 =%10.3f cost1 =%10.3f\n", it1->Cost[0], it1->Cost[1]);
        int len = it1->sizeOfAddO_origin;
        int t=0;
        for (int j=0; j<len; j+=3){
            fprintf(outFile, "%4d ", t);
            switch (t % 4){
                case 0: fprintf(outFile, " N  "); break;
                case 1: fprintf(outFile, " CA "); break;
                case 2: fprintf(outFile, " C  "); break;
                case 3: fprintf(outFile, " O  "); break;
            }
            fprintf(outFile, "%8.3f %8.3f %8.3f\n",
                    it1->addO_origin[j], it1->addO_origin[j+1], it1->addO_origin[j+2]);
            t++;
        }
        k++;
        fclose(outFile);
    }
}

void printTimeForDebug(int choice, time_t start, time_t finish){
    const char *timeFile = catStrStr(draftAddress, "time.txt");
    ofstream outfile(timeFile, ios::app);
    long duration;
    duration = finish-start;
    switch (choice){
        case 0: outfile << "TOTAL Time: \n"; break;
        case 1: outfile << "Cost 1 : "; break;
        case 2: outfile << "Cost 2 : "; break;
        case 3: outfile << "Transfer Time: "; break;
        case 4: outfile << "Input Time: "; break;
        case 5: outfile << "Adaption Time: "; break;
        case 6: outfile << "IntoRep Time: "; break;
        default: break;
    }
    int seconds = duration % 60;
    duration /= 60;
    int minutes = duration % 60;
    duration /= 60;
    int hours = duration;

    outfile << hours << ":" << minutes << ":" << seconds << endl;
    outfile.close();
}

void printRep(myRep & rep, int len){
    int repNum=0;
    repNum++;
    char *filename;
    filename = catStrIntStr("/home/ws/zzZyj/data/rep/rep", repNum, ".txt");
    ofstream outFile(filename);
    myRep ::iterator it2;
    int rr=0;
    for (it2=rep.begin(); it2!= rep.end(); it2++){
        rr++;
        outFile <<  rr << " cost[0]=" << it2->Cost[0] << " cost[1]=" << it2->Cost[1] << endl;
        outFile << "Following is the x,y,z coordinary: "<< endl;
        for (int i=0; i<len; i+=3){
            outFile<< i/3+1 << " " << it2->addO_origin[i] << " " << it2->addO_origin[i+1] << " " << it2->addO_origin[i+2]<< endl;
        }
        outFile << "------------------------------------------------------------------"<< endl;
    }
    delete filename;
    outFile.close();
}

// file stream
/*void printAnswer(myRep & rep){
    myRep ::iterator it1;
    int k=0;
    for (it1=rep.begin(); it1!= rep.end(); it1++){
        char *filename;
        filename = catStrIntStr("/home/ws/zzZyj/data/answer/rep", k, ".txt");
        ofstream outFile(filename);
        outFile << " cost0=" << it1->Cost[0] << " cost1=" << it1->Cost[1] << endl;
        int len = it1->sizeOfAddO_origin;
        int t=0;
        for (int j=0; j<len; j+=3){
            outFile << t << " ";
            switch (t % 4){
                case 0: outFile << "N "; break;
                case 1: outFile << "CA "; break;
                case 2: outFile << "C "; break;
                case 3: outFile << "O "; break;
            }
            outFile << it1->addO_origin[j] << " " << it1->addO_origin[j+1] << " " << it1->addO_origin[j+2] << endl;
            t++;
        }
        outFile << endl;
        k++;
    }
}*/
