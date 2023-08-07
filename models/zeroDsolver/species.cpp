#include<species.H>

namespace mflo_species
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real molwts[NUM_SPECIES]={one};
    AMREX_GPU_DEVICE_MANAGED amrex::Real advect_flags[NUM_SPECIES];
    AMREX_GPU_DEVICE_MANAGED amrex::Real rateconsts[NUM_REACTIONS];
    AMREX_GPU_DEVICE_MANAGED amrex::Real gamma_spec[NUM_GAS_SPECIES+1];
    AMREX_GPU_DEVICE_MANAGED amrex::Real rxnarray[NUM_REACTIONS][NUM_SPECIES];

    void init()
    {
        specnames[CH4_ID] = "CH4";
        specnames[Ar_ID] = "Ar";
        specnames[C2H6_ID] = "C2H6";
        specnames[C2H4_ID] = "C2H4";
        specnames[C2H2_ID] = "C2H2";
        specnames[C3H8_ID] = "C3H8";
        specnames[C3H6_ID] = "C3H6";
        specnames[O2_ID] = "O2";
        specnames[H2O_ID] = "H2O";
        specnames[H2O2_ID] = "H2O2";
        specnames[CH2O_ID] = "CH2O";
        specnames[CO_ID] = "CO";
        specnames[CO2_ID] = "CO2";
        specnames[H2_ID] = "H2";
        specnames[HO2_ID] = "HO2";
        specnames[O_ID] = "O";
        specnames[OH_ID] = "OH";
        specnames[CH3_ID] = "CH3";
        specnames[C2H5_ID] = "C2H5";
        specnames[C2H3_ID] = "C2H3";
        specnames[C3H7_ID] = "C3H7";
        specnames[CH3O_ID] = "CH3O";
        specnames[CHO_ID] = "CHO";
        specnames[H_ID] = "H";
        specnames[S1_ID] = "S1";
        specnames[OADS1_ID] = "OADS1";
        specnames[OHADS1_ID] = "OHADS1";
        specnames[H2OADS1_ID] = "H2OADS1";
        specnames[CH3OADS1_ID] = "CH3OADS1";
        specnames[CH2OADS1_ID] = "CH2OADS1";
        specnames[CHOADS1_ID] = "CHOADS1";
        specnames[COADS1_ID] = "COADS1";
        specnames[CO2ADS1_ID] = "CO2ADS1";

        //kg/mol
        molwts[CH4_ID] = 0.01604;
        molwts[Ar_ID] = 0.03995;
        molwts[C2H6_ID] = 0.03007;
        molwts[C2H4_ID] = 0.02805;
        molwts[C2H2_ID] = 0.02604;
        molwts[C3H8_ID] = 0.04410;
        molwts[C3H6_ID] = 0.04208;
        molwts[O2_ID] = 0.03200;
        molwts[H2O_ID] = 0.01801;
        molwts[H2O2_ID] = 0.03401;
        molwts[CH2O_ID] = 0.03003;
        molwts[CO_ID] = 0.02801;
        molwts[CO2_ID] = 0.04401;
        molwts[H2_ID] = 0.00202;
        molwts[HO2_ID] = 0.03301;
        molwts[O_ID] = 0.01600;
        molwts[OH_ID] = 0.01701;
        molwts[CH3_ID] = 0.01503;
        molwts[C2H5_ID] = 0.02906;
        molwts[C2H3_ID] = 0.02705;
        molwts[C3H7_ID] = 0.04309;
        molwts[CH3O_ID] = 0.03103;
        molwts[CHO_ID] = 0.02902;
        molwts[H_ID] = 0.00101;
        molwts[S1_ID] = 0.00000;
        molwts[OADS1_ID] = 0.01600;
        molwts[OHADS1_ID] = 0.01701;
        molwts[H2OADS1_ID] = 0.01801;
        molwts[CH3OADS1_ID] = 0.03103;
        molwts[CH2OADS1_ID] = 0.03003;
        molwts[CHOADS1_ID] = 0.02902;
        molwts[COADS1_ID] = 0.02801;
        molwts[CO2ADS1_ID] = 0.04401;

        //gamma for each species
        gamma_spec[CH4_ID] = 1.13;
        gamma_spec[Ar_ID] = 1.60;
        gamma_spec[C2H6_ID] = 1.07;
        gamma_spec[C2H4_ID] = 1.10;
        gamma_spec[C2H2_ID] = 1.14;
        gamma_spec[C3H8_ID] = 1.05;
        gamma_spec[C3H6_ID] = 1.06;
        gamma_spec[O2_ID] = 1.31;
        gamma_spec[H2O_ID] = 1.25;
        gamma_spec[H2O2_ID] = 1.15;
        gamma_spec[CH2O_ID] = 1.15;
        gamma_spec[CO_ID] = 1.33;
        gamma_spec[CO2_ID] = 1.06;
        gamma_spec[H2_ID] = 1.00;
        gamma_spec[HO2_ID] = 1.21;
        gamma_spec[O_ID] = 1.66;
        gamma_spec[OH_ID] = 1.37;
        gamma_spec[CH3_ID] = 1.16;
        gamma_spec[C2H5_ID] = 1.08;
        gamma_spec[C2H3_ID] = 1.12;
        gamma_spec[C3H7_ID] = 1.05;
        gamma_spec[CH3O_ID] = 1.12;
        gamma_spec[CHO_ID] = 1.21;
        gamma_spec[H_ID] = 1.67;
        //background gas
        gamma_spec[NUM_GAS_SPECIES]  = GAMMA_BG_GAS;



        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            advect_flags[sp]=one;
        }
        advect_flags[S1_ID] = zeroval;
        advect_flags[OADS1_ID] = zeroval;
        advect_flags[OHADS1_ID] = zeroval;
        advect_flags[H2OADS1_ID] = zeroval;
        advect_flags[CH3OADS1_ID] = zeroval;
        advect_flags[CH2OADS1_ID] = zeroval;
        advect_flags[CHOADS1_ID] = zeroval;
        advect_flags[COADS1_ID] = zeroval;
        advect_flags[CO2ADS1_ID] = zeroval;
        for(int i=0;i<NUM_REACTIONS;i++)
        {
            for(int j=0;j<NUM_SPECIES;j++)
            {
                rxnarray[i][j]=zeroval;
            }
        }
        //Reaction0
        rxnarray[0][CH4_ID] = -1;
        rxnarray[0][O2_ID] = -1;
        rxnarray[0][HO2_ID] = 1;
        rxnarray[0][CH3_ID] = 1;
        rateconsts[0] = 1.24e-03;

        //Reaction1
        rxnarray[1][CH4_ID] = -1;
        rxnarray[1][H2_ID] = 1;
        rxnarray[1][CH3_ID] = 1;
        rxnarray[1][H_ID] = -1;
        rateconsts[1] = 5.71e+05;

        //Reaction2
        rxnarray[2][CH4_ID] = -1;
        rxnarray[2][O_ID] = -1;
        rxnarray[2][OH_ID] = 1;
        rxnarray[2][CH3_ID] = 1;
        rateconsts[2] = 2.38e+07;

        //Reaction3
        rxnarray[3][CH4_ID] = -1;
        rxnarray[3][H2O_ID] = 1;
        rxnarray[3][OH_ID] = -1;
        rxnarray[3][CH3_ID] = 1;
        rateconsts[3] = 5.69e+06;

        //Reaction4
        rxnarray[4][CH4_ID] = -1;
        rxnarray[4][H2O2_ID] = 1;
        rxnarray[4][HO2_ID] = -1;
        rxnarray[4][CH3_ID] = 1;
        rateconsts[4] = 3.29e+02;

        //Reaction5
        rxnarray[5][O2_ID] = -1;
        rxnarray[5][O_ID] = 1;
        rxnarray[5][CH3_ID] = -1;
        rxnarray[5][CH3O_ID] = 1;
        rateconsts[5] = 1.94e+01;

        //Reaction6
        rxnarray[6][O2_ID] = -1;
        rxnarray[6][CH2O_ID] = 1;
        rxnarray[6][OH_ID] = 1;
        rxnarray[6][CH3_ID] = -1;
        rateconsts[6] = 2.34e+02;

        //Reaction7
        rxnarray[7][HO2_ID] = -1;
        rxnarray[7][OH_ID] = 1;
        rxnarray[7][CH3_ID] = -1;
        rxnarray[7][CH3O_ID] = 1;
        rateconsts[7] = 8.85e+07;

        //Reaction8
        rxnarray[8][C2H6_ID] = 1;
        rxnarray[8][CH3_ID] = -2;
        rateconsts[8] = 6.5e+07;

        //Reaction9
        rxnarray[9][CH2O_ID] = 1;
        rxnarray[9][CH3O_ID] = -1;
        rxnarray[9][H_ID] = 1;
        rateconsts[9] = 3.46e+08;

        //Reaction10
        rxnarray[10][H2O_ID] = 1;
        rxnarray[10][CH2O_ID] = -1;
        rxnarray[10][OH_ID] = -1;
        rxnarray[10][CHO_ID] = 1;
        rateconsts[10] = 3.22e+08;

        //Reaction11
        rxnarray[11][H2O2_ID] = 1;
        rxnarray[11][CH2O_ID] = -1;
        rxnarray[11][HO2_ID] = -1;
        rxnarray[11][CHO_ID] = 1;
        rateconsts[11] = 3.73e+04;

        //Reaction12
        rxnarray[12][CH4_ID] = 1;
        rxnarray[12][CH2O_ID] = -1;
        rxnarray[12][CH3_ID] = -1;
        rxnarray[12][CHO_ID] = 1;
        rateconsts[12] = 3.69E+06;

        //Reaction13
        rxnarray[13][CO_ID] = 1;
        rxnarray[13][CHO_ID] = -1;
        rxnarray[13][H_ID] = 1;
        rateconsts[13] = 1.45E+06;

        //Reaction14
        rxnarray[14][O2_ID] = -1;
        rxnarray[14][CO_ID] = 1;
        rxnarray[14][HO2_ID] = 1;
        rxnarray[14][CHO_ID] = -1;
        rateconsts[14] = 1.71E+05;

        //Reaction15
        rxnarray[15][CO_ID] = -1;
        rxnarray[15][CO2_ID] = 1;
        rxnarray[15][HO2_ID] = -1;
        rxnarray[15][OH_ID] = 1;
        rateconsts[15] = 1.02E+03;

        //Reaction16
        rxnarray[16][C2H6_ID] = -1;
        rxnarray[16][H2_ID] = 1;
        rxnarray[16][C2H5_ID] = 1;
        rxnarray[16][H_ID] = -1;
        rateconsts[16] = 2.09E+06;

        //Reaction17
        rxnarray[17][C2H6_ID] = -1;
        rxnarray[17][H2O_ID] = 1;
        rxnarray[17][OH_ID] = -1;
        rxnarray[17][C2H5_ID] = 1;
        rateconsts[17] = 7.25E+07;

        //Reaction18
        rxnarray[18][CH4_ID] = 1;
        rxnarray[18][C2H6_ID] = -1;
        rxnarray[18][CH3_ID] = -1;
        rxnarray[18][C2H5_ID] = 1;
        rateconsts[18] = 1.18E+04;

        //Reaction19
        rxnarray[19][CH2O_ID] = 1;
        rxnarray[19][HO2_ID] = -1;
        rxnarray[19][OH_ID] = 1;
        rxnarray[19][CH3_ID] = 1;
        rxnarray[19][C2H5_ID] = -1;
        rateconsts[19] = 9.48E+06;

        //Reaction20
        rxnarray[20][C2H4_ID] = 1;
        rxnarray[20][C2H5_ID] = -1;
        rxnarray[20][H_ID] = 1;
        rateconsts[20] = 1.64E+05;

        //Reaction21
        rxnarray[21][C2H4_ID] = 1;
        rxnarray[21][O2_ID] = -1;
        rxnarray[21][HO2_ID] = 1;
        rxnarray[21][C2H5_ID] = -1;
        rateconsts[21] = 1.22E+04;

        //Reaction22
        rxnarray[22][C2H4_ID] = -1;
        rxnarray[22][O2_ID] = -1;
        rxnarray[22][HO2_ID] = 1;
        rxnarray[22][C2H3_ID] = 1;
        rateconsts[22] = 1.17E-01;

        //Reaction23
        rxnarray[23][C2H4_ID] = -1;
        rxnarray[23][H2_ID] = 1;
        rxnarray[23][C2H3_ID] = 1;
        rxnarray[23][H_ID] = -1;
        rateconsts[23] = 9.90E+05;

        //Reaction24
        rxnarray[24][C2H4_ID] = -1;
        rxnarray[24][H2O_ID] = 1;
        rxnarray[24][OH_ID] = -1;
        rxnarray[24][C2H3_ID] = 1;
        rateconsts[24] = 3.35E+06;

        //Reaction25
        rxnarray[25][CH4_ID] = 1;
        rxnarray[25][C2H4_ID] = -1;
        rxnarray[25][CH3_ID] = -1;
        rxnarray[25][C2H3_ID] = 1;
        rateconsts[25] = 4.69E+02;

        //Reaction26
        rxnarray[26][C2H4_ID] = -1;
        rxnarray[26][CH2O_ID] = 1;
        rxnarray[26][OH_ID] = -1;
        rxnarray[26][CH3_ID] = 1;
        rateconsts[26] = 2.72E+06;

        //Reaction27
        rxnarray[27][C2H2_ID] = 1;
        rxnarray[27][C2H3_ID] = -1;
        rxnarray[27][H_ID] = 1;
        rateconsts[27] = 1.18E+06;

        //Reaction28
        rxnarray[28][C2H2_ID] = 1;
        rxnarray[28][O2_ID] = -1;
        rxnarray[28][HO2_ID] = 1;
        rxnarray[28][C2H3_ID] = -1;
        rateconsts[28] = 5.00E+06;

        //Reaction29
        rxnarray[29][O2_ID] = -1;
        rxnarray[29][CH2O_ID] = 1;
        rxnarray[29][C2H3_ID] = -1;
        rxnarray[29][CHO_ID] = 1;
        rateconsts[29] = 5.50E+06;

        //Reaction30
        rxnarray[30][C3H8_ID] = 1;
        rxnarray[30][CH3_ID] = -1;
        rxnarray[30][C2H5_ID] = -1;
        rateconsts[30] = 8.00E+06;

        //Reaction31
        rxnarray[31][C3H8_ID] = -1;
        rxnarray[31][H2_ID] = 1;
        rxnarray[31][C3H7_ID] = 1;
        rxnarray[31][H_ID] = -1;
        rateconsts[31] = 2.09E+07;

        //Reaction32
        rxnarray[32][C2H4_ID] = -1;
        rxnarray[32][CH3_ID] = -1;
        rxnarray[32][C3H7_ID] = 1;
        rateconsts[32] = 9.92E+03;

        //Reaction33
        rxnarray[33][C3H6_ID] = 1;
        rxnarray[33][C3H7_ID] = -1;
        rxnarray[33][H_ID] = 1;
        rateconsts[33] = 1.62E+07;

        //Reaction34
        rxnarray[34][O2_ID] = -1;
        rxnarray[34][O_ID] = 1;
        rxnarray[34][OH_ID] = 1;
        rxnarray[34][H_ID] = -1;
        rateconsts[34] = 5.66E+04;

        //Reaction35
        rxnarray[35][O2_ID] = -1;
        rxnarray[35][HO2_ID] = 1;
        rxnarray[35][H_ID] = -1;
        rateconsts[35] = 1.39E+05;

        //Reaction36
        rxnarray[36][O2_ID] = 1;
        rxnarray[36][HO2_ID] = -2;
        rxnarray[36][OH_ID] = 2;
        rateconsts[36] = 2.00E+06;

        //Reaction37
        rxnarray[37][H2O2_ID] = -1;
        rxnarray[37][OH_ID] = 2;
        rateconsts[37] = 8.40E+00;

        //Reaction38
        rxnarray[38][C2H6_ID] = -1;
        rxnarray[38][C2H5_ID] = 1;
        rxnarray[38][H_ID] = 1;
        rateconsts[38] = 1.88E-03;

        //Reaction39
        rxnarray[39][CH4_ID] = 1;
        rxnarray[39][O2_ID] = 1;
        rxnarray[39][HO2_ID] = -1;
        rxnarray[39][CH3_ID] = -1;
        rateconsts[39] = 1.07E+07;

        //Reaction40
        rxnarray[40][CH4_ID] = 1;
        rxnarray[40][H2_ID] = -1;
        rxnarray[40][CH3_ID] = -1;
        rxnarray[40][H_ID] = 1;
        rateconsts[40] = 0.00000; //actually zero

        //Reaction41
        rxnarray[41][CH4_ID] = 1;
        rxnarray[41][O_ID] = 1;
        rxnarray[41][OH_ID] = -1;
        rxnarray[41][CH3_ID] = -1;
        rateconsts[41] = 1.28E+06;

        //Reaction42
        rxnarray[42][CH4_ID] = 1;
        rxnarray[42][H2O_ID] = -1;
        rxnarray[42][OH_ID] = 1;
        rxnarray[42][CH3_ID] = -1;
        rateconsts[42] = 7.46E+02;

        //Reaction43
        rxnarray[43][CH4_ID] = 1;
        rxnarray[43][H2O2_ID] = -1;
        rxnarray[43][HO2_ID] = 1;
        rxnarray[43][CH3_ID] = -1;
        rateconsts[43] = 2.48E+05;

        //Reaction44
        rxnarray[44][O2_ID] = 1;
        rxnarray[44][O_ID] = -1;
        rxnarray[44][CH3_ID] = 1;
        rxnarray[44][CH3O_ID] = -1;
        rateconsts[44] = 1.03E+08;

        //Reaction45
        rxnarray[45][O2_ID] = 1;
        rxnarray[45][CH2O_ID] = -1;
        rxnarray[45][OH_ID] = -1;
        rxnarray[45][CH3_ID] = 1;
        rateconsts[45] = 8.19e-10;

        //Reaction46
        rxnarray[46][HO2_ID] = 1;
        rxnarray[46][OH_ID] = -1;
        rxnarray[46][CH3_ID] = 1;
        rxnarray[46][CH3O_ID] = -1;
        rateconsts[46] = 2.92E+03;

        //Reaction47
        rxnarray[47][C2H6_ID] = -1;
        rxnarray[47][CH3_ID] = 2;
        rateconsts[47] = 2.24E-03;

        //Reaction48
        rxnarray[48][CH2O_ID] = -1;
        rxnarray[48][CH3O_ID] = 1;
        rxnarray[48][H_ID] = -1;
        rateconsts[48] = 1.21E+07;

        //Reaction49
        rxnarray[49][H2O_ID] = -1;
        rxnarray[49][CH2O_ID] = 1;
        rxnarray[49][OH_ID] = 1;
        rxnarray[49][CHO_ID] = -1;
        rateconsts[49] = 5.13E+01;

        //Reaction50
        rxnarray[50][H2O2_ID] = -1;
        rxnarray[50][CH2O_ID] = 1;
        rxnarray[50][HO2_ID] = 1;
        rxnarray[50][CHO_ID] = -1;
        rateconsts[50] = 3.41E+04;

        //Reaction51
        rxnarray[51][CH4_ID] = -1;
        rxnarray[51][CH2O_ID] = 1;
        rxnarray[51][CH3_ID] = 1;
        rxnarray[51][CHO_ID] = -1;
        rateconsts[51] = 4.48E+03;

        //Reaction52
        rxnarray[52][CO_ID] = -1;
        rxnarray[52][CHO_ID] = 1;
        rxnarray[52][H_ID] = -1;
        rateconsts[52] = 3.27E+04;

        //Reaction53
        rxnarray[53][O2_ID] = 1;
        rxnarray[53][CO_ID] = -1;
        rxnarray[53][HO2_ID] = -1;
        rxnarray[53][CHO_ID] = 1;
        rateconsts[53] = 1.18E-02;

        //Reaction54
        rxnarray[54][CO_ID] = 1;
        rxnarray[54][CO2_ID] = -1;
        rxnarray[54][HO2_ID] = 1;
        rxnarray[54][OH_ID] = -1;
        rateconsts[54] = 0.00000; //actually zero

        //Reaction55
        rxnarray[55][C2H6_ID] = 1;
        rxnarray[55][H2_ID] = -1;
        rxnarray[55][C2H5_ID] = -1;
        rxnarray[55][H_ID] = 1;
        rateconsts[55] = 0.00000; //actually zero

        //Reaction56
        rxnarray[56][C2H6_ID] = 1;
        rxnarray[56][H2O_ID] = -1;
        rxnarray[56][OH_ID] = 1;
        rxnarray[56][C2H5_ID] = -1;
        rateconsts[56] = 2.87E+02;

        //Reaction57
        rxnarray[57][CH4_ID] = -1;
        rxnarray[57][C2H6_ID] = 1;
        rxnarray[57][CH3_ID] = 1;
        rxnarray[57][C2H5_ID] = -1;
        rateconsts[57] = 3.57E+02;

        //Reaction58
        rxnarray[58][CH2O_ID] = -1;
        rxnarray[58][HO2_ID] = 1;
        rxnarray[58][OH_ID] = -1;
        rxnarray[58][CH3_ID] = -1;
        rxnarray[58][C2H5_ID] = 1;
        rateconsts[58] = 3.70E-03;

        //Reaction59
        rxnarray[59][C2H4_ID] = -1;
        rxnarray[59][C2H5_ID] = 1;
        rxnarray[59][H_ID] = -1;
        rateconsts[59] = 2.08E+08;

        //Reaction60
        rxnarray[60][C2H4_ID] = -1;
        rxnarray[60][O2_ID] = 1;
        rxnarray[60][HO2_ID] = -1;
        rxnarray[60][C2H5_ID] = 1;
        rateconsts[60] = 4.74E+01;

        //Reaction61
        rxnarray[61][C2H4_ID] = 1;
        rxnarray[61][O2_ID] = 1;
        rxnarray[61][HO2_ID] = -1;
        rxnarray[61][C2H3_ID] = -1;
        rateconsts[61] = 4.34E+09;

        //Reaction62
        rxnarray[62][C2H4_ID] = 1;
        rxnarray[62][H2_ID] = -1;
        rxnarray[62][C2H3_ID] = -1;
        rxnarray[62][H_ID] = 1;
        rateconsts[62] = 0.00000; //actually zero

        //Reaction63
        rxnarray[63][C2H4_ID] = 1;
        rxnarray[63][H2O_ID] = -1;
        rxnarray[63][OH_ID] = 1;
        rxnarray[63][C2H3_ID] = -1;
        rateconsts[63] = 1.89E+03;

        //Reaction64
        rxnarray[64][CH4_ID] = -1;
        rxnarray[64][C2H4_ID] = 1;
        rxnarray[64][CH3_ID] = 1;
        rxnarray[64][C2H3_ID] = -1;
        rateconsts[64] = 6.23E+07;

        //Reaction65
        rxnarray[65][C2H4_ID] = 1;
        rxnarray[65][CH2O_ID] = -1;
        rxnarray[65][OH_ID] = 1;
        rxnarray[65][CH3_ID] = -1;
        rateconsts[65] = 6.65E+02;

        //Reaction66
        rxnarray[66][C2H2_ID] = -1;
        rxnarray[66][C2H3_ID] = 1;
        rxnarray[66][H_ID] = -1;
        rateconsts[66] = 1.53E+09;

        //Reaction67
        rxnarray[67][C2H2_ID] = -1;
        rxnarray[67][O2_ID] = 1;
        rxnarray[67][HO2_ID] = -1;
        rxnarray[67][C2H3_ID] = 1;
        rateconsts[67] = 1.97E+04;

        //Reaction68
        rxnarray[68][O2_ID] = 1;
        rxnarray[68][CH2O_ID] = -1;
        rxnarray[68][C2H3_ID] = 1;
        rxnarray[68][CHO_ID] = -1;
        rateconsts[68] = 1.33e-12;

        //Reaction69
        rxnarray[69][C3H8_ID] = -1;
        rxnarray[69][CH3_ID] = 1;
        rxnarray[69][C2H5_ID] = 1;
        rateconsts[69] = 2.52E-03;

        //Reaction70
        rxnarray[70][C3H8_ID] = 1;
        rxnarray[70][H2_ID] = -1;
        rxnarray[70][C3H7_ID] = -1;
        rxnarray[70][H_ID] = 1;
        rateconsts[70] = 0.00000; //actually zero

        //Reaction71
        rxnarray[71][C2H4_ID] = 1;
        rxnarray[71][CH3_ID] = 1;
        rxnarray[71][C3H7_ID] = -1;
        rateconsts[71] = 11.91E+06;

        //Reaction72
        rxnarray[72][C3H6_ID] = -1;
        rxnarray[72][C3H7_ID] = 1;
        rxnarray[72][H_ID] = -1;
        rateconsts[72] = 8.58E+08;

        //Reaction73
        rxnarray[73][O2_ID] = 1;
        rxnarray[73][O_ID] = -1;
        rxnarray[73][OH_ID] = -1;
        rxnarray[73][H_ID] = 1;
        rateconsts[73] = 1.14E+07;

        //Reaction74
        rxnarray[74][O2_ID] = 1;
        rxnarray[74][HO2_ID] = -1;
        rxnarray[74][H_ID] = 1;
        rateconsts[74] = 4.24E-01;

        //Reaction75
        rxnarray[75][O2_ID] = -1;
        rxnarray[75][HO2_ID] = 2;
        rxnarray[75][OH_ID] = -2;
        rateconsts[75] = 8.22E+02;

        //Reaction76
        rxnarray[76][H2O2_ID] = 1;
        rxnarray[76][OH_ID] = -2;
        rateconsts[76] = 1.25E+09;

        //Reaction77
        rxnarray[77][C2H6_ID] = 1;
        rxnarray[77][C2H5_ID] = -1;
        rxnarray[77][H_ID] = -1;
        rateconsts[77] = 1.61e11;
/*
        //Reaction78
        rxnarray[78][O2_ID] = -1;
        rxnarray[78][S1_ID] = -2;
        rxnarray[78][OADS1_ID] = 2;
        rateconsts[78] = 7.85e16;

        //Reaction79
        rxnarray[79][CH4_ID] = -1;
        rxnarray[79][CH3_ID] = 1;
        rxnarray[79][OADS1_ID] = -1;
        rxnarray[79][OHADS1_ID] = 1;
        rateconsts[79] = 6.19E-01;

        //Reaction80
        rxnarray[80][C2H4_ID] = -1;
        rxnarray[80][C2H3_ID] = 1;
        rxnarray[80][OADS1_ID] = -1;
        rxnarray[80][OHADS1_ID] = 1;
        rateconsts[80] = 2.73E+00;

        //Reaction81
        rxnarray[81][C2H6_ID] = -1;
        rxnarray[81][C2H5_ID] = 1;
        rxnarray[81][OADS1_ID] = -1;
        rxnarray[81][OHADS1_ID] = 1;
        rateconsts[81] = 4.40E-02;

        //Reaction82
        rxnarray[82][OADS1_ID] = 1;
        rxnarray[82][OHADS1_ID] = -2;
        rxnarray[82][H2OADS1_ID] = 1;
        rateconsts[82] = 4.32E+05;

        //Reaction83
        rxnarray[83][H2O_ID] = 1;
        rxnarray[83][S1_ID] = 1;
        rxnarray[83][H2OADS1_ID] = -1;
        rateconsts[83] = 3.59E+10;

        //Reaction84
        rxnarray[84][CH3_ID] = -1;
        rxnarray[84][OADS1_ID] = -1;
        rxnarray[84][CH3OADS1_ID] = 1;
        rateconsts[84] = 5.30E+07;

        //Reaction85
        rxnarray[85][OADS1_ID] = -1;
        rxnarray[85][OHADS1_ID] = 1;
        rxnarray[85][CH3OADS1_ID] = -1;
        rxnarray[85][CH2OADS1_ID] = 1;
        rateconsts[85] = 2.25E+27;

        //Reaction86
        rxnarray[86][OADS1_ID] = -1;
        rxnarray[86][OHADS1_ID] = 1;
        rxnarray[86][CH2OADS1_ID] = -1;
        rxnarray[86][CHOADS1_ID] = 1;
        rateconsts[86] = 2.73E+13;

        //Reaction87
        rxnarray[87][OADS1_ID] = -1;
        rxnarray[87][OHADS1_ID] = 1;
        rxnarray[87][CHOADS1_ID] = -1;
        rxnarray[87][COADS1_ID] = 1;
        rateconsts[87] = 3.11E+14;

        //Reaction88
        rxnarray[88][S1_ID] = 1;
        rxnarray[88][OADS1_ID] = -1;
        rxnarray[88][COADS1_ID] = -1;
        rxnarray[88][CO2ADS1_ID] = 1;
        rateconsts[88] = 2.44E+27;

        //Reaction89
        rxnarray[89][CO_ID] = -1;
        rxnarray[89][S1_ID] = -1;
        rxnarray[89][COADS1_ID] = 1;
        rateconsts[89] = 4.57E+08;

        //Reaction90
        rxnarray[90][CO2_ID] = -1;
        rxnarray[90][S1_ID] = -1;
        rxnarray[90][CO2ADS1_ID] = 1;
        rateconsts[90] = 1.35E+10;

        //Reaction91
        rxnarray[91][O2_ID] = 3;
        rxnarray[91][H2O_ID] = 2;
        rxnarray[91][HO2_ID] = -4;
        rateconsts[91] = 1.00E-02;

        //Reaction92
        rxnarray[92][O2_ID] = 1;
        rxnarray[92][S1_ID] = 2;
        rxnarray[92][OADS1_ID] = -2;
        rateconsts[92] = 5.32E+09;

        //Reaction93
        rxnarray[93][CH4_ID] = 1;
        rxnarray[93][CH3_ID] = -1;
        rxnarray[93][OADS1_ID] = 1;
        rxnarray[93][OHADS1_ID] = -1;
        rateconsts[93] = 8.72E+02;

        //Reaction94
        rxnarray[94][C2H4_ID] = 1;
        rxnarray[94][C2H3_ID] = -1;
        rxnarray[94][OADS1_ID] = 1;
        rxnarray[94][OHADS1_ID] = -1;
        rateconsts[94] = 3.56E+02;

        //Reaction95
        rxnarray[95][C2H6_ID] = 1;
        rxnarray[95][C2H5_ID] = -1;
        rxnarray[95][OADS1_ID] = 1;
        rxnarray[95][OHADS1_ID] = -1;
        rateconsts[95] = 1.34E+03;

        //Reaction96
        rxnarray[96][OADS1_ID] = -1;
        rxnarray[96][OHADS1_ID] = 2;
        rxnarray[96][H2OADS1_ID] = -1;
        rateconsts[96] = 2.10E+10;

        //Reaction97
        rxnarray[97][H2O_ID] = -1;
        rxnarray[97][S1_ID] = -1;
        rxnarray[97][H2OADS1_ID] = 1;
        rateconsts[97] = 1.10E+12;

        //Reaction98
        rxnarray[98][CH3_ID] = 1;
        rxnarray[98][OADS1_ID] = 1;
        rxnarray[98][CH3OADS1_ID] = -1;
        rateconsts[98] = 7.25E+00;

        //Reaction99
        rxnarray[99][OADS1_ID] = 1;
        rxnarray[99][OHADS1_ID] = -1;
        rxnarray[99][CH3OADS1_ID] = 1;
        rxnarray[99][CH2OADS1_ID] = -1;
        rateconsts[99] = 1.85E+07;

        //Reaction100
        rxnarray[100][OADS1_ID] = 1;
        rxnarray[100][OHADS1_ID] = -1;
        rxnarray[100][CH2OADS1_ID] = 1;
        rxnarray[100][CHOADS1_ID] = -1;
        rateconsts[100] = 3.11E+09;

        //Reaction101
        rxnarray[101][OADS1_ID] = 1;
        rxnarray[101][OHADS1_ID] = -1;
        rxnarray[101][CHOADS1_ID] = 1;
        rxnarray[101][COADS1_ID] = -1;
        rateconsts[101] = 2.63E+08;

        //Reaction102
        rxnarray[102][S1_ID] = -1;
        rxnarray[102][OADS1_ID] = 1;
        rxnarray[102][COADS1_ID] = 1;
        rxnarray[102][CO2ADS1_ID] = -1;
        rateconsts[102] = 4.57E+04;

        //Reaction103
        rxnarray[103][CO_ID] = 1;
        rxnarray[103][S1_ID] = 1;
        rxnarray[103][COADS1_ID] = -1;
        rateconsts[103] = 1.20E+09;

        //Reaction104
        rxnarray[104][CO2_ID] = 1;
        rxnarray[104][S1_ID] = 1;
        rxnarray[104][CO2ADS1_ID] = -1;
        rateconsts[104] = 6.51E+04;

        */

    }    
    void close()
    {
        specnames.clear();
    }
    int find_id(std::string specname)
    {
        int loc=-1;
        auto it=std::find(specnames.begin(),specnames.end(),specname);
        if(it != specnames.end())
        {
            loc=it-specnames.begin();
        }
        return(loc);
    }
}