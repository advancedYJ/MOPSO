#include "MOPSO.h"

char *catStrIntStr(const char *c1, int ss,const char *c2){
    int len1 = strlen(c1), len2 = strlen(c2);
    char *str3;
    str3 = new char [len1+len2+10];
    strcpy(str3, c1);
    char tmp[100];
    sprintf(tmp, "%d", ss);
    strcat(str3, tmp);
    strcat(str3, c2);

    return str3;
}

char *catStrIntStr(const char *c1, const char *c2, int ss, const char *c3){
    const char *str1 = catStrStr(c1, c2);
    return catStrIntStr(str1, ss, c3);
}

char *catStrStr(const char *c1, const char *c2){
    int len1 = strlen(c1), len2 = strlen(c2);
    char *str3;
    str3 = new char [len1+len2+10];
    strcpy(str3, c1);
    strcat(str3,c2);
    return str3;
}

char *catStrStr(const char *c1, const char *c2, const char *c3){
    char *str = catStrStr(c1, c2);
    return catStrStr(str, c3);
}

char *catStrStr(const char *c1, const char *c2, const char *c3, const char *c4){
    const char *str = catStrStr(c1, c2);
    return catStrStr(str, c3, c4);
}

char *catStrStr(const char *c1, const char *c2, const char *c3, const char *c4, const char *c5){
    const char *str = catStrStr(c1, c2);
    return catStrStr(str, c3, c4, c5);
}

char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6){
    const char *str = catStrStr(c1, c2);
    return catStrStr(str, c3, c4, c5, c6);
}

char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6, CONCH *c7){
    const char *str = catStrStr(c1, c2);
    return catStrStr(str, c3, c4, c5, c6, c7);
}

void cpyDoubleArray(double *& p1, double *& p2, int n){                //  p1,p2 is Position,  copy p2 to p1
    for (int i=0; i<n; i++)
        p1[i]  = p2[i];
}

void checkV(double *& tmp, int n){
    for (int i=0; i<n; i++)
        if (tmp[i] > VelMax)        tmp[i] = VelMax;
        else if (tmp[i] < -VelMax)      tmp[i] = -VelMax;
}

void checkP(double *& tmp, double *& tmpV, int n){
    for (int i=0; i<n; i++)
        if (tmp[i] > VarMax) {  tmp[i] = VarMax;  tmpV[i] = -tmpV[i]; }
        else if (tmp[i] < VarMin)   {    tmp[i] = VarMin;    tmpV[i] = -tmpV[i]; }
    //  if the particle reach the boudary, v need to inverse
}

double getInertiaWeight(double it, double MaxIt){
    return wMax- (wMax-wMin)*(it+1)/MaxIt;
}

int rouletteWheel(myRep & rep){
    int sz = rep.size();
    return rand() % sz;
    //return 0;
}

double dis1(double x, double y, double z){
    double ans = (x-1)*(x-1) + (y-5)*(y-5) + (z+12)*(z+12);
    return sqrt(ans);
}

double dis2(double x, double y, double z){
    double ans = (x+1)*(x+1) + (y-3)*(y-3) + (z+7)*(z+7);
    return sqrt(ans);
}

char * getRelativeAA(char abbreviation){
    char *str;
    str = new char [3];
    switch (abbreviation){
        case 'A':   strcpy(str, "ALA"); break;
        case 'R':   strcpy(str, "ARG"); break;
        case 'N':   strcpy(str, "ASN"); break;
        case 'D':   strcpy(str, "ASP"); break;
        case 'C':   strcpy(str, "CYS"); break;
        case 'Q':   strcpy(str, "GLN"); break;
        case 'E':   strcpy(str, "GLU"); break;
        case 'G':   strcpy(str, "GLY"); break;
        case 'H':   strcpy(str, "HIS"); break;
        case 'I':   strcpy(str, "ILE"); break;
        case 'L':   strcpy(str, "LEU"); break;
        case 'K':   strcpy(str, "LYS"); break;
        case 'M':   strcpy(str, "MET"); break;
        case 'F':   strcpy(str, "PHE"); break;
        case 'P':   strcpy(str, "PRO"); break;
        case 'S':   strcpy(str, "SER"); break;
        case 'T':   strcpy(str, "THR"); break;
        case 'W':   strcpy(str, "TRP"); break;
        case 'Y':   strcpy(str, "TYR"); break;
        default :   strcpy(str, "VAL"); break;
    }
    return str;
}

double getCost0(int i) {
    //remove("/home/ws/zzZyj/data/temp1.pdb"); remove("default.sc");
    /*system("/home/ws/GL/rosetta_source/bin/score.linuxgccrelease -database \
                   /home/ws/GL/rosetta_database/ \
                   -s /home/ws/zzZyj/data/energyFile/tempFile/temp_1.pdb \
                   -out:file:scorefile /home/ws/zzZyj/data/default0.sc");*/
    const char *temp_File = catStrIntStr(energyFileAddress, "tempFile/temp_", i, ".pdb");
    const char *defaultFile = catStrIntStr(energyFileAddress, "/defaultFile/default", i, ".sc");
    /*const char *command = catStrStr("/home/ws/GL/rosetta_source/bin/score.linuxgccrelease \
                                            -database /home/ws/GL/rosetta_database/ -s ",
                    temp_File, " -out:file:scorefile ", defaultFile);*/
    const char *command = catStrStr(scoreAddress, " -database ", databaseAddress," -s ",
                                    temp_File, " -out:file:scorefile ", defaultFile);
    system(command);
    ifstream infile;
    openFile(defaultFile, infile);
    char buffer[10000], strAns[10000];
    double ans=0;
    infile.getline(buffer, 10000);
    infile >> buffer >> strAns;
    int len = strlen(strAns);
    if (strAns[0]=='n' || (len>=2 && strAns[1]=='n'))
        ans = INF;
    else
        ans = atof(strAns);
    infile.close();
    removeFile(defaultFile);
    return ans;
}
/*
double getCost1(int i){
    //remove("/home/ws/GL/mybin/temp1.pdb");  remove("QUACKout.txt");
    //system("cp /home/ws/zzZyj/data/temp1.pdb /home/ws/GL/mybin/temp1.pdb");
    const char *command01 = catStrIntStr(energyFileAddress, "tempFile/temp_", i, ".pdb");
    const char *command02 = catStrIntStr(mybinAddress, "temp_", i, ".pdb");
    const char *command0 = catStrStr("cp ", command01, " ", command02);
    system(command0);

    double ans=0;
    while (ans==0){
        int lineNum=0;
        char buffer[10000];
        //system("/home/ws/GL/mybin/calcquarkenergy /home/ws/GL/mybin/ /home/ws/GL/mybin/temp1.pdb>>QUACKout.txt");
        const char *str = catStrStr(mybinAddress, "calcquarkenergy ", mybinAddress, " ",
                                    mybinAddress, "temp_");
        const char *command1 = catStrIntStr(str, i, ".pdb>>");
        const char *command2 = catStrIntStr(energyFileAddress, "QUACKoutFile/QUACKout", i, ".txt");
        const char *command = catStrStr(command1, command2);
        system(command);
        ifstream infile;
        openFile(command2, infile);
        while (infile.getline(buffer, 10000)) lineNum++;
        infile.close();
        openFile(command2, infile);
        for (int i=0; i<lineNum-1; i++) infile.getline(buffer, 10000);
        infile >> ans;
        infile.close();
        removeFile(command2);
        delete command;  delete command1; delete command2;
    }
    removeFile(command01);
    removeFile(command02);

    delete command0;  delete command01; delete command02;

    return ans;
}
 */