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
        specnames[Ar_ID] = "Ar";
        specnames[CH4_ID] = "CH4";
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

        //kg/mol
        molwts[Ar_ID] = 0.03995;
        molwts[CH4_ID] = 0.01604;
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
        molwts[O_ID] = 15.99900;
        molwts[OH_ID] = 0.01701;
        molwts[CH3_ID] = 0.01503;
        molwts[C2H5_ID] = 0.02906;
        molwts[C2H3_ID] = 0.02705;
        molwts[C3H7_ID] = 0.04309;
        molwts[CH3O_ID] = 0.03103;
        molwts[CHO_ID] = 0.02902;
        molwts[H_ID] = 0.00101;
        molwts[S1_ID] = 0.00000;

        //gamma for each species
        gamma_spec[Ar_ID] = 1.6666666669999999;
        gamma_spec[CH4_ID] = 1.1292790129999999;
        gamma_spec[C2H6_ID] = 1.0719076220000001;
        gamma_spec[C2H4_ID] = 1.0962566470000001;
        gamma_spec[C2H2_ID] = 1.1381456559999998;
        gamma_spec[C3H8_ID] = 1.049710951;
        gamma_spec[C3H6_ID] = 1.060734871;
        gamma_spec[O2_ID] = 1.311803071;
        gamma_spec[H2O_ID] = 1.250600705;
        gamma_spec[H2O2_ID] = 1.149667264;
        gamma_spec[CH2O_ID] = 1.1532652559999998;
        gamma_spec[CO_ID] = 1.33252325;
        gamma_spec[CO2_ID] = 1.057985308;
        gamma_spec[H2_ID] = 1.001386814;
        gamma_spec[HO2_ID] = 1.210657844;
        gamma_spec[O_ID] = 1.6599627030000002;
        gamma_spec[OH_ID] = 1.3706028719999999;
        gamma_spec[CH3_ID] = 1.1628463740000001;
        gamma_spec[C2H5_ID] = 1.084061698;
        gamma_spec[C2H3_ID] = 1.118041235;
        gamma_spec[C3H7_ID] = 1.054438897;
        gamma_spec[CH3O_ID] = 1.118955465;
        gamma_spec[CHO_ID] = 1.207677183;
        gamma_spec[H_ID] = 1.6666666669999999;
        //background gas
        gamma_spec[NUM_GAS_SPECIES]  = GAMMA_BG_GAS;



        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            advect_flags[sp]=one;
        }
        advect_flags[S1_ID] = zeroval;
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
        rateconsts[0] = 0.00124;

        //Reaction1
        rxnarray[1][CH4_ID] = -1;
        rxnarray[1][H2_ID] = 1;
        rxnarray[1][CH3_ID] = 1;
        rxnarray[1][H_ID] = -1;
        rateconsts[1] = 571151.53230;

        //Reaction2
        rxnarray[2][CH4_ID] = -1;
        rxnarray[2][O_ID] = -1;
        rxnarray[2][OH_ID] = 1;
        rxnarray[2][CH3_ID] = 1;
        rateconsts[2] = 23802698.03000;

        //Reaction3
        rxnarray[3][CH4_ID] = -1;
        rxnarray[3][H2O_ID] = 1;
        rxnarray[3][OH_ID] = 1;
        rxnarray[3][CH3_ID] = 1;
        rateconsts[3] = 5698968.10500;

        //Reaction4
        rxnarray[4][CH4_ID] = -1;
        rxnarray[4][H2O2_ID] = 1;
        rxnarray[4][HO2_ID] = -1;
        rxnarray[4][CH3_ID] = 1;
        rateconsts[4] = 329.30206;

        //Reaction5
        rxnarray[5][O2_ID] = -1;
        rxnarray[5][O_ID] = 1;
        rxnarray[5][CH3_ID] = -1;
        rxnarray[5][CH3O_ID] = 1;
        rateconsts[5] = 19.49173;

        //Reaction6
        rxnarray[6][O2_ID] = -1;
        rxnarray[6][CH2O_ID] = 1;
        rxnarray[6][OH_ID] = 1;
        rxnarray[6][CH3_ID] = -1;
        rateconsts[6] = 234.14861;

        //Reaction7
        rxnarray[7][HO2_ID] = -1;
        rxnarray[7][OH_ID] = 1;
        rxnarray[7][CH3_ID] = -1;
        rxnarray[7][CH3O_ID] = 1;
        rateconsts[7] = 88500000.00000;

        //Reaction8
        rxnarray[8][C2H6_ID] = 1;
        rxnarray[8][CH3_ID] = -2;
        rateconsts[8] = 65000000.00000;

        //Reaction9
        rxnarray[9][CH2O_ID] = 1;
        rxnarray[9][CH3O_ID] = -1;
        rxnarray[9][H_ID] = 1;
        rateconsts[9] = 34700784400.00000;

        //Reaction10
        rxnarray[10][H2O_ID] = 1;
        rxnarray[10][CH2O_ID] = -1;
        rxnarray[10][OH_ID] = -1;
        rxnarray[10][CHO_ID] = 1;
        rateconsts[10] = 322221809.40000;

        //Reaction11
        rxnarray[11][CH4_ID] = 1;
        rxnarray[11][CH2O_ID] = -1;
        rxnarray[11][HO2_ID] = -1;
        rxnarray[11][CHO_ID] = 1;
        rateconsts[11] = 37309.97944;

        //Reaction12
        rxnarray[12][CH4_ID] = 1;
        rxnarray[12][CH2O_ID] = -1;
        rxnarray[12][CH3_ID] = -1;
        rxnarray[12][CHO_ID] = 1;
        rateconsts[12] = 3691484.88300;

        //Reaction13
        rxnarray[13][CO_ID] = 1;
        rxnarray[13][CHO_ID] = -1;
        rxnarray[13][H_ID] = 1;
        rateconsts[13] = 1449716.84100;

        //Reaction14
        rxnarray[14][O2_ID] = -1;
        rxnarray[14][CO_ID] = 1;
        rxnarray[14][HO2_ID] = 1;
        rxnarray[14][CHO_ID] = -1;
        rateconsts[14] = 171000.00000;

        //Reaction15
        rxnarray[15][CO_ID] = -1;
        rxnarray[15][CO2_ID] = 1;
        rxnarray[15][HO2_ID] = -1;
        rxnarray[15][OH_ID] = 1;
        rateconsts[15] = 1019.40940;

        //Reaction16
        rxnarray[16][C2H6_ID] = -1;
        rxnarray[16][H2_ID] = 1;
        rxnarray[16][C2H5_ID] = 1;
        rxnarray[16][H_ID] = -1;
        rateconsts[16] = 2086978.31600;

        //Reaction17
        rxnarray[17][C2H6_ID] = -1;
        rxnarray[17][H2O_ID] = 1;
        rxnarray[17][OH_ID] = -1;
        rxnarray[17][C2H5_ID] = 1;
        rateconsts[17] = 72493547.94000;

        //Reaction18
        rxnarray[18][CH4_ID] = 1;
        rxnarray[18][C2H6_ID] = -1;
        rxnarray[18][CH3_ID] = -1;
        rxnarray[18][C2H5_ID] = 1;
        rateconsts[18] = 11847.66722;

        //Reaction19
        rxnarray[19][CH2O_ID] = 1;
        rxnarray[19][HO2_ID] = -1;
        rxnarray[19][OH_ID] = 1;
        rxnarray[19][CH3_ID] = 1;
        rxnarray[19][C2H5_ID] = -1;
        rateconsts[19] = 9480000.00000;

        //Reaction20
        rxnarray[20][C2H4_ID] = 1;
        rxnarray[20][C2H5_ID] = -1;
        rxnarray[20][H_ID] = 1;
        rateconsts[20] = 164221.48020;

        //Reaction21
        rxnarray[21][C2H4_ID] = 1;
        rxnarray[21][O2_ID] = -1;
        rxnarray[21][HO2_ID] = 1;
        rxnarray[21][C2H5_ID] = -1;
        rateconsts[21] = 12208.66989;

        //Reaction22
        rxnarray[22][C2H4_ID] = -1;
        rxnarray[22][O2_ID] = -1;
        rxnarray[22][HO2_ID] = 1;
        rxnarray[22][C2H3_ID] = 1;
        rateconsts[22] = 0.11716;

        //Reaction23
        rxnarray[23][C2H4_ID] = -1;
        rxnarray[23][H2_ID] = 1;
        rxnarray[23][C2H3_ID] = 1;
        rxnarray[23][H_ID] = -1;
        rateconsts[23] = 990967.27760;

        //Reaction24
        rxnarray[24][C2H4_ID] = -1;
        rxnarray[24][H2O_ID] = 1;
        rxnarray[24][OH_ID] = -1;
        rxnarray[24][C2H3_ID] = 1;
        rateconsts[24] = 3355077.05100;

        //Reaction25
        rxnarray[25][CH4_ID] = 1;
        rxnarray[25][C2H4_ID] = -1;
        rxnarray[25][CH3_ID] = -1;
        rxnarray[25][C2H3_ID] = 1;
        rateconsts[25] = 469.44285;

        //Reaction26
        rxnarray[26][C2H4_ID] = -1;
        rxnarray[26][CH2O_ID] = 1;
        rxnarray[26][OH_ID] = -1;
        rxnarray[26][CH3_ID] = 1;
        rateconsts[26] = 2720000.00000;

        //Reaction27
        rxnarray[27][C2H2_ID] = 1;
        rxnarray[27][C2H3_ID] = -1;
        rxnarray[27][H_ID] = 1;
        rateconsts[27] = 1187707.71300;

        //Reaction28
        rxnarray[28][C2H2_ID] = 1;
        rxnarray[28][O2_ID] = -1;
        rxnarray[28][HO2_ID] = 1;
        rxnarray[28][C2H3_ID] = -1;
        rateconsts[28] = 5000000.00000;

        //Reaction29
        rxnarray[29][O2_ID] = -1;
        rxnarray[29][CH2O_ID] = 1;
        rxnarray[29][C2H3_ID] = -1;
        rxnarray[29][CHO_ID] = 1;
        rateconsts[29] = 5500000.00000;

        //Reaction30
        rxnarray[30][C3H8_ID] = 1;
        rxnarray[30][CH3_ID] = -1;
        rxnarray[30][C2H5_ID] = -1;
        rateconsts[30] = 8000000.00000;

        //Reaction31
        rxnarray[31][C3H8_ID] = -1;
        rxnarray[31][H2_ID] = 1;
        rxnarray[31][C3H7_ID] = 1;
        rxnarray[31][H_ID] = -1;
        rateconsts[31] = 20916784.21000;

        //Reaction32
        rxnarray[32][C2H4_ID] = -1;
        rxnarray[32][CH3_ID] = -1;
        rxnarray[32][C3H7_ID] = 1;
        rateconsts[32] = 9920.58993;

        //Reaction33
        rxnarray[33][C3H6_ID] = 1;
        rxnarray[33][C3H7_ID] = -1;
        rxnarray[33][H_ID] = 1;
        rateconsts[33] = 16276901.71000;

        //Reaction34
        rxnarray[34][O2_ID] = -1;
        rxnarray[34][O_ID] = 1;
        rxnarray[34][OH_ID] = 1;
        rxnarray[34][H_ID] = -1;
        rateconsts[34] = 56660.89581;

        //Reaction35
        rxnarray[35][O2_ID] = -1;
        rxnarray[35][HO2_ID] = 1;
        rxnarray[35][H_ID] = -1;
        rateconsts[35] = 139000.00000;

        //Reaction36
        rxnarray[36][O2_ID] = 1;
        rxnarray[36][HO2_ID] = -2;
        rxnarray[36][OH_ID] = 2;
        rateconsts[36] = 2000000.00000;

        //Reaction37
        rxnarray[37][H2O2_ID] = -1;
        rxnarray[37][OH_ID] = 2;
        rateconsts[37] = 8.42473;

        //Reaction38
        rxnarray[38][C2H6_ID] = -1;
        rxnarray[38][C2H5_ID] = 1;
        rxnarray[38][H_ID] = 1;
        rateconsts[38] = 0.00189;

        //Reaction39
        rxnarray[39][CH4_ID] = 1;
        rxnarray[39][O2_ID] = 1;
        rxnarray[39][HO2_ID] = -1;
        rxnarray[39][CH3_ID] = -1;
        rateconsts[39] = 10724388.76000;

        //Reaction40
        rxnarray[40][CH4_ID] = 1;
        rxnarray[40][H2_ID] = -1;
        rxnarray[40][CH3_ID] = -1;
        rxnarray[40][H_ID] = 1;
        rateconsts[40] = 0.00000;

        //Reaction41
        rxnarray[41][CH4_ID] = 1;
        rxnarray[41][O_ID] = 1;
        rxnarray[41][OH_ID] = -1;
        rxnarray[41][CH3_ID] = -1;
        rateconsts[41] = 1284416.27800;

        //Reaction42
        rxnarray[42][CH4_ID] = 1;
        rxnarray[42][H2O_ID] = -1;
        rxnarray[42][OH_ID] = -1;
        rxnarray[42][CH3_ID] = -1;
        rateconsts[42] = 747.63243;

        //Reaction43
        rxnarray[43][CH4_ID] = 1;
        rxnarray[43][H2O2_ID] = -1;
        rxnarray[43][HO2_ID] = 1;
        rxnarray[43][CH3_ID] = -1;
        rateconsts[43] = 247692.18620;

        //Reaction44
        rxnarray[44][O2_ID] = 1;
        rxnarray[44][O_ID] = -1;
        rxnarray[44][CH3_ID] = 1;
        rxnarray[44][CH3O_ID] = -1;
        rateconsts[44] = 102729755.60000;

        //Reaction45
        rxnarray[45][O2_ID] = 1;
        rxnarray[45][CH2O_ID] = -1;
        rxnarray[45][OH_ID] = -1;
        rxnarray[45][CH3_ID] = 1;
        rateconsts[45] = 0.00000;

        //Reaction46
        rxnarray[46][HO2_ID] = 1;
        rxnarray[46][OH_ID] = -1;
        rxnarray[46][CH3_ID] = 1;
        rxnarray[46][CH3O_ID] = -1;
        rateconsts[46] = 2921.47777;

        //Reaction47
        rxnarray[47][C2H6_ID] = -1;
        rxnarray[47][CH3_ID] = 2;
        rateconsts[47] = 0.00225;

        //Reaction48
        rxnarray[48][CH2O_ID] = -1;
        rxnarray[48][CH3O_ID] = 1;
        rxnarray[48][H_ID] = -1;
        rateconsts[48] = 1207707236.00000;

        //Reaction49
        rxnarray[49][H2O_ID] = -1;
        rxnarray[49][CH2O_ID] = 1;
        rxnarray[49][OH_ID] = 1;
        rxnarray[49][CHO_ID] = -1;
        rateconsts[49] = 51.42969;

        //Reaction50
        rxnarray[50][CH4_ID] = -1;
        rxnarray[50][CH2O_ID] = 1;
        rxnarray[50][HO2_ID] = 1;
        rxnarray[50][CHO_ID] = -1;
        rateconsts[50] = 34143.64925;

        //Reaction51
        rxnarray[51][CH4_ID] = -1;
        rxnarray[51][CH2O_ID] = 1;
        rxnarray[51][CH3_ID] = 1;
        rxnarray[51][CHO_ID] = -1;
        rateconsts[51] = 4491.25958;

        //Reaction52
        rxnarray[52][CO_ID] = -1;
        rxnarray[52][CHO_ID] = 1;
        rxnarray[52][H_ID] = -1;
        rateconsts[52] = 32655.85650;

        //Reaction53
        rxnarray[53][O2_ID] = 1;
        rxnarray[53][CO_ID] = -1;
        rxnarray[53][HO2_ID] = -1;
        rxnarray[53][CHO_ID] = 1;
        rateconsts[53] = 0.01180;

        //Reaction54
        rxnarray[54][CO_ID] = 1;
        rxnarray[54][CO2_ID] = -1;
        rxnarray[54][HO2_ID] = 1;
        rxnarray[54][OH_ID] = -1;
        rateconsts[54] = 0.00000;

        //Reaction55
        rxnarray[55][C2H6_ID] = 1;
        rxnarray[55][H2_ID] = -1;
        rxnarray[55][C2H5_ID] = -1;
        rxnarray[55][H_ID] = 1;
        rateconsts[55] = 0.00000;

        //Reaction56
        rxnarray[56][C2H6_ID] = 1;
        rxnarray[56][H2O_ID] = -1;
        rxnarray[56][OH_ID] = 1;
        rxnarray[56][C2H5_ID] = -1;
        rateconsts[56] = 287.21102;

        //Reaction57
        rxnarray[57][CH4_ID] = -1;
        rxnarray[57][C2H6_ID] = 1;
        rxnarray[57][CH3_ID] = 1;
        rxnarray[57][C2H5_ID] = -1;
        rateconsts[57] = 357.80194;

        //Reaction58
        rxnarray[58][CH2O_ID] = -1;
        rxnarray[58][HO2_ID] = 1;
        rxnarray[58][OH_ID] = -1;
        rxnarray[58][CH3_ID] = -1;
        rxnarray[58][C2H5_ID] = 1;
        rateconsts[58] = 0.00370;

        //Reaction59
        rxnarray[59][C2H4_ID] = -1;
        rxnarray[59][C2H5_ID] = 1;
        rxnarray[59][H_ID] = -1;
        rateconsts[59] = 208334841.40000;

        //Reaction60
        rxnarray[60][C2H4_ID] = -1;
        rxnarray[60][O2_ID] = 1;
        rxnarray[60][HO2_ID] = -1;
        rxnarray[60][C2H5_ID] = 1;
        rateconsts[60] = 47.44488;

        //Reaction61
        rxnarray[61][C2H4_ID] = 1;
        rxnarray[61][O2_ID] = 1;
        rxnarray[61][HO2_ID] = -1;
        rxnarray[61][C2H3_ID] = -1;
        rateconsts[61] = 4336801408.00000;

        //Reaction62
        rxnarray[62][C2H4_ID] = 1;
        rxnarray[62][H2_ID] = -1;
        rxnarray[62][C2H3_ID] = -1;
        rxnarray[62][H_ID] = 1;
        rateconsts[62] = 0.00000;

        //Reaction63
        rxnarray[63][C2H4_ID] = 1;
        rxnarray[63][H2O_ID] = -1;
        rxnarray[63][OH_ID] = 1;
        rxnarray[63][C2H3_ID] = -1;
        rateconsts[63] = 1891.19216;

        //Reaction64
        rxnarray[64][CH4_ID] = -1;
        rxnarray[64][C2H4_ID] = 1;
        rxnarray[64][CH3_ID] = 1;
        rxnarray[64][C2H3_ID] = -1;
        rateconsts[64] = 62263090.40000;

        //Reaction65
        rxnarray[65][C2H4_ID] = 1;
        rxnarray[65][CH2O_ID] = -1;
        rxnarray[65][OH_ID] = 1;
        rxnarray[65][CH3_ID] = -1;
        rateconsts[65] = 665.45423;

        //Reaction66
        rxnarray[66][C2H2_ID] = -1;
        rxnarray[66][C2H3_ID] = 1;
        rxnarray[66][H_ID] = -1;
        rateconsts[66] = 1526702637.00000;

        //Reaction67
        rxnarray[67][C2H2_ID] = -1;
        rxnarray[67][O2_ID] = 1;
        rxnarray[67][HO2_ID] = -1;
        rxnarray[67][C2H3_ID] = 1;
        rateconsts[67] = 19688.10698;

        //Reaction68
        rxnarray[68][O2_ID] = 1;
        rxnarray[68][CH2O_ID] = -1;
        rxnarray[68][C2H3_ID] = 1;
        rxnarray[68][CHO_ID] = -1;
        rateconsts[68] = 0.00000;

        //Reaction69
        rxnarray[69][C3H8_ID] = -1;
        rxnarray[69][CH3_ID] = 1;
        rxnarray[69][C2H5_ID] = 1;
        rateconsts[69] = 0.00254;

        //Reaction70
        rxnarray[70][C3H8_ID] = 1;
        rxnarray[70][H2_ID] = -1;
        rxnarray[70][C3H7_ID] = -1;
        rxnarray[70][H_ID] = 1;
        rateconsts[70] = 0.00000;

        //Reaction71
        rxnarray[71][C2H4_ID] = 1;
        rxnarray[71][CH3_ID] = 1;
        rxnarray[71][C3H7_ID] = -1;
        rateconsts[71] = 1909954.60300;

        //Reaction72
        rxnarray[72][C3H6_ID] = -1;
        rxnarray[72][C3H7_ID] = 1;
        rxnarray[72][H_ID] = -1;
        rateconsts[72] = 857943173.50000;

        //Reaction73
        rxnarray[73][O2_ID] = 1;
        rxnarray[73][O_ID] = -1;
        rxnarray[73][OH_ID] = -1;
        rxnarray[73][H_ID] = 1;
        rateconsts[73] = 11377229.45000;

        //Reaction74
        rxnarray[74][O2_ID] = 1;
        rxnarray[74][HO2_ID] = -1;
        rxnarray[74][H_ID] = 1;
        rateconsts[74] = 0.42580;

        //Reaction75
        rxnarray[75][O2_ID] = -1;
        rxnarray[75][HO2_ID] = 2;
        rxnarray[75][OH_ID] = -2;
        rateconsts[75] = 821.12153;

        //Reaction76
        rxnarray[76][H2O2_ID] = 1;
        rxnarray[76][OH_ID] = -2;
        rateconsts[76] = 1257277733.00000;

        //Reaction77
        rxnarray[77][C2H6_ID] = 1;
        rxnarray[77][C2H5_ID] = -1;
        rxnarray[77][H_ID] = -1;
        rateconsts[77] = 161000000000.00000;

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