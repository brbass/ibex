#ifndef Linear_Algebra_hh
#define Linear_Algebra_hh

#include "Check.hh"

#include <vector>

namespace Linear_Algebra
{
    template<class T> void linear_solve_1(std::vector<T> const &a,
                                          std::vector<T> const &b,
                                          std::vector<T> &x)
    {
        x[0]=b[0]/a[0];
        return;
    }
    
    template<class T> void symmetric_linear_solve_1(std::vector<T> const &a,
                                                    std::vector<T> const &b,
                                                    std::vector<T> &x)
    {
        x[0]=b[0]/a[0];
        return;
    }
    
    template<class T> void linear_solve_2(std::vector<T> const &a,
                                          std::vector<T> const &b,
                                          std::vector<T> &x)
    {
        T o7=a[3];
        T o3=a[1];
        T o4=a[2];
        T o6=a[0];
        T o11=b[0];
        T o13=b[1];
        x[0]=(-(o13*o3)+o11*o7)/(-(o3*o4)+o6*o7);
        x[1]=(o11*o4-o13*o6)/(o3*o4-o6*o7);
        return;
    }

    template<class T> void symmetric_linear_solve_2(std::vector<T> const &a,
                                                    std::vector<T> const &b,
                                                    std::vector<T> &x)
    {
        T o7=a[2];
        T o3=a[1];
        T o4=o3*o3;
        T o6=a[0];
        T o11=b[0];
        T o13=b[1];
        x[0]=(-(o13*o3)+o11*o7)/(-o4+o6*o7);
        x[1]=(o11*o3-o13*o6)/(o4-o6*o7);
        return;
    }

    template<class T> void linear_solve_3(std::vector<T> const &a,
                                          std::vector<T> const &b,
                                          std::vector<T> &x)
    {
        T o5=a[6];
        T o3=a[2];
        T o8=a[5];
        T o11=a[7];
        T o7=a[1];
        T o10=a[3];
        T o13=a[0];
        T o4=a[4];
        T o15=a[8];
        T o20=b[0];
        T o23=b[1];
        T o26=b[2];
        T o6=o3*o4*o5;
        T o9=-(o5*o7*o8);
        T o12=-(o10*o11*o3);
        T o14=o11*o13*o8;
        T o16=o10*o15*o7;
        T o17=-(o13*o15*o4);
        T o18=o12+o14+o16+o17+o6+o9;
        T o19=1/o18;
        x[0]=o19*(-(o11*o23*o3)-o15*o20*o4+o26*o3*o4+o15*o23*o7+o11*o20*o8-o26*o7*o8);
        x[1]=o19*(o10*o15*o20-o13*o15*o23-o10*o26*o3+o23*o3*o5+o13*o26*o8-o20*o5*o8);
        x[2]=o19*(-(o10*o11*o20)+o11*o13*o23-o13*o26*o4+o20*o4*o5+o10*o26*o7-o23*o5*o7);
        return;
    }

    template<class T> void symmetric_linear_solve_3(std::vector<T> const &a,
                                                    std::vector<T> const &b,
                                                    std::vector<T> &x)
    {
        T o3=a[2];
        T o4=a[3];
        T o8=a[4];
        T o10=a[1];
        T o12=a[0];
        T o15=a[5];
        T o13=o8*o8;
        T o20=b[0];
        T o23=b[1];
        T o26=b[2];
        T o5=o4*o4;
        T o6=o3*o5;
        T o7=o3*o3;
        T o9=-(o7*o8);
        T o11=-(o10*o4*o8);
        T o14=o12*o13;
        T o16=o10*o15*o3;
        T o17=-(o12*o15*o4);
        T o18=o11+o14+o16+o17+o6+o9;
        T o19=1/o18;
        x[0]=o19*(o13*o20+o10*o15*o23-o15*o20*o4+o26*o3*o4-o10*o26*o8-o23*o3*o8);
        x[1]=o19*(-(o12*o15*o23)+o15*o20*o3+o23*o3*o4-o26*o7+o12*o26*o8-o20*o4*o8);
        x[2]=(-(o10*o26*o3)+o10*o23*o4+o12*o26*o4-o20*o5-o12*o23*o8+o20*o3*o8)/(-(o12*o13)-o10*o15*o3+o12*o15*o4-o3*o5+o10*o4*o8+o7*o8);        
    }

    template<class T> void linear_solve_4(std::vector<T> const &a,
                                          std::vector<T> const &b,
                                          std::vector<T> &x)
    {
        T o5=a[9];
        T o6=a[12];
        T o3=a[3];
        T o9=a[7];
        T o12=a[10];
        T o8=a[2];
        T o11=a[5];
        T o14=a[1];
        T o4=a[6];
        T o16=a[11];
        T o19=a[8];
        T o20=a[13];
        T o23=a[4];
        T o25=a[0];
        T o29=a[14];
        T o36=a[15];
        T o45=b[0];
        T o52=b[1];
        T o59=b[2];
        T o66=b[3];
        T o7=-(o3*o4*o5*o6);
        T o10=o5*o6*o8*o9;
        T o13=o11*o12*o3*o6;
        T o15=-(o12*o14*o6*o9);
        T o17=-(o11*o16*o6*o8);
        T o18=o14*o16*o4*o6;
        T o21=o19*o20*o3*o4;
        T o22=-(o19*o20*o8*o9);
        T o24=-(o12*o20*o23*o3);
        T o26=o12*o20*o25*o9;
        T o27=o16*o20*o23*o8;
        T o28=-(o16*o20*o25*o4);
        T o30=-(o11*o19*o29*o3);
        T o31=o14*o19*o29*o9;
        T o32=o23*o29*o3*o5;
        T o33=-(o25*o29*o5*o9);
        T o34=-(o14*o16*o23*o29);
        T o35=o11*o16*o25*o29;
        T o37=o11*o19*o36*o8;
        T o38=-(o14*o19*o36*o4);
        T o39=-(o23*o36*o5*o8);
        T o40=o25*o36*o4*o5;
        T o41=o12*o14*o23*o36;
        T o42=-(o11*o12*o25*o36);
        T o43=o10+o13+o15+o17+o18+o21+o22+o24+o26+o27+o28+o30+o31+o32+o33+o34+o35+o37+o38+o39+o40+o41+o42+o7;
        T o44=1/o43;
        T o134=-(o14*o23);
        T o135=o11*o25;
        T o136=o134+o135;
        T o142=-(o23*o8);
        T o143=o25*o4;
        T o144=o142+o143;
        T o130=-(o14*o19);
        T o131=o25*o5;
        T o132=o130+o131;
        T o127=-(o23*o3);
        T o128=o25*o9;
        T o129=o127+o128;
        T o145=-(o14*o6);
        T o146=o20*o25;
        T o147=o145+o146;
        T o148=-(o144*o147);
        T o149=-(o6*o8);
        T o150=o25*o29;
        T o151=o149+o150;
        T o152=o136*o151;
        T o153=o148+o152;
        T o155=-(o132*o144);
        T o156=-(o19*o8);
        T o157=o12*o25;
        T o158=o156+o157;
        T o159=o136*o158;
        T o160=o155+o159;
        T o170=-(o23*o45);
        T o171=o25*o52;
        T o172=o170+o171;
        x[0]=o44*(o11*o16*o29*o45-o11*o12*o36*o45-o16*o20*o4*o45+o36*o4*o45*o5-o14*o16*o29*o52-o12*o20*o3*o52+o12*o14*o36*o52+o29*o3*o5*o52-o11*o29*o3*o59+o20*o3*o4*o59-o14*o36*o4*o59+o11*o12*o3*o66+o14*o16*o4*o66-o3*o4*o5*o66+o16*o20*o52*o8-o36*o5*o52*o8+o11*o36*o59*o8-o11*o16*o66*o8+o12*o20*o45*o9-o29*o45*o5*o9+o14*o29*o59*o9-o12*o14*o66*o9-o20*o59*o8*o9+o5*o66*o8*o9);
        x[1]=o44*(-(o16*o23*o29*o45)+o12*o23*o36*o45-o19*o36*o4*o45+o16*o25*o29*o52-o19*o29*o3*o52-o12*o25*o36*o52+o23*o29*o3*o59+o25*o36*o4*o59+o16*o4*o45*o6+o12*o3*o52*o6-o3*o4*o59*o6-o12*o23*o3*o66-o16*o25*o4*o66+o19*o3*o4*o66+o19*o36*o52*o8-o23*o36*o59*o8-o16*o52*o6*o8+o16*o23*o66*o8+o19*o29*o45*o9-o25*o29*o59*o9-o12*o45*o6*o9+o12*o25*o66*o9+o59*o6*o8*o9-o19*o66*o8*o9);
        x[2]=o44*(o16*o20*o23*o45+o11*o19*o36*o45-o23*o36*o45*o5-o16*o20*o25*o52+o19*o20*o3*o52-o14*o19*o36*o52+o25*o36*o5*o52-o20*o23*o3*o59+o14*o23*o36*o59-o11*o25*o36*o59-o11*o16*o45*o6+o14*o16*o52*o6-o3*o5*o52*o6+o11*o3*o59*o6-o14*o16*o23*o66+o11*o16*o25*o66-o11*o19*o3*o66+o23*o3*o5*o66-o19*o20*o45*o9+o20*o25*o59*o9+o45*o5*o6*o9-o14*o59*o6*o9+o14*o19*o66*o9-o25*o5*o66*o9);
        x[3]=(-(o153*(-(o132*o172)+o136*(-(o19*o45)+o25*o59)))+o160*(-(o147*o172)+o136*(-(o45*o6)+o25*o66)))/(-(o153*(-(o129*o132)+o136*(o16*o25-o19*o3)))+o160*(-(o129*o147)+o136*(o25*o36-o3*o6)));
        return;
    }

    template<class T> void symmetric_linear_solve_4(std::vector<T> const &a,
                                                    std::vector<T> const &b,
                                                    std::vector<T> &x)
    {
        T o5=a[6];
        T o3=a[3];
        T o4=a[5];
        T o12=a[7];
        T o8=a[2];
        T o6=o5*o5;
        T o20=o12*o12;
        T o14=a[4];
        T o19=o3*o3;
        T o24=a[8];
        T o17=a[1];
        T o22=a[0];
        T o32=o24*o24;
        T o11=o4*o4;
        T o35=a[9];
        T o44=b[0];
        T o51=b[1];
        T o58=b[2];
        T o65=b[3];
        T o7=-(o3*o4*o6);
        T o9=pow(o5,3.);
        T o10=o8*o9;
        T o13=o11*o12*o3;
        T o15=o12*o14*o3*o5;
        T o16=-(o12*o4*o5*o8);
        T o18=-(o12*o17*o6);
        T o21=-(o19*o20);
        T o23=o20*o22*o5;
        T o25=-(o14*o24*o3*o4);
        T o26=o19*o24*o5;
        T o27=-(o14*o24*o5*o8);
        T o28=2.*o17*o24*o4*o5;
        T o29=-(o22*o24*o6);
        T o30=o12*o24*o3*o8;
        T o31=-(o12*o22*o24*o4);
        T o33=-(o17*o3*o32);
        T o34=o14*o22*o32;
        T o36=o14*o35*o4*o8;
        T o37=-(o11*o17*o35);
        T o38=-(o3*o35*o5*o8);
        T o39=o22*o35*o4*o5;
        T o40=o12*o17*o3*o35;
        T o41=-(o12*o14*o22*o35);
        T o42=o10+o13+o15+o16+o18+o21+o23+o25+o26+o27+o28+o29+o30+o31+o33+o34+o36+o37+o38+o39+o40+o41+o7;
        T o43=1/o42;
        T o126=o22*o5;
        T o131=-(o17*o3);
        T o132=o14*o22;
        T o133=o131+o132;
        T o135=o22*o24;
        T o139=-(o3*o8);
        T o140=o22*o4;
        T o141=o139+o140;
        T o128=-(o17*o4);
        T o129=o126+o128;
        T o143=o12*o22;
        T o125=-o19;
        T o127=o125+o126;
        T o142=-(o17*o5);
        T o144=o142+o143;
        T o145=-(o141*o144);
        T o146=-(o5*o8);
        T o147=o135+o146;
        T o148=o133*o147;
        T o149=o145+o148;
        T o151=-(o129*o141);
        T o152=-(o4*o8);
        T o153=o143+o152;
        T o154=o133*o153;
        T o155=o151+o154;
        T o165=-(o3*o44);
        T o166=o22*o51;
        T o167=o165+o166;
        x[0]=o43*(o14*o32*o44-o12*o14*o35*o44-o12*o24*o4*o44+o20*o44*o5+o35*o4*o44*o5-o20*o3*o51-o17*o32*o51+o12*o17*o35*o51+o24*o3*o5*o51-o14*o24*o3*o58+o12*o3*o4*o58-o17*o35*o4*o58+o17*o24*o5*o58-o24*o44*o6+o12*o14*o3*o65+o17*o24*o4*o65-o12*o17*o5*o65-o3*o4*o5*o65+o12*o24*o51*o8-o35*o5*o51*o8+o14*o35*o58*o8-o12*o5*o58*o8-o14*o24*o65*o8+o6*o65*o8);
        x[1]=o43*(-(o3*o32*o44)-o11*o35*o44+o12*o3*o35*o44+2.*o24*o4*o44*o5+o22*o32*o51-o12*o22*o35*o51-o24*o3*o4*o51+o12*o3*o5*o51+o19*o24*o58+o22*o35*o4*o58-o22*o24*o5*o58-o3*o4*o5*o58-o12*o44*o6-o12*o19*o65+o11*o3*o65-o22*o24*o4*o65+o12*o22*o5*o65+o35*o4*o51*o8-o24*o5*o51*o8-o3*o35*o58*o8+o58*o6*o8+o24*o3*o65*o8-o4*o5*o65*o8);
        x[2]=o43*(o12*o24*o3*o44+o14*o35*o4*o44-o14*o24*o44*o5-o3*o35*o44*o5-o12*o4*o44*o5-o12*o22*o24*o51+o12*o3*o4*o51-o17*o35*o4*o51+o17*o24*o5*o51+o22*o35*o5*o51-o12*o19*o58-o14*o22*o35*o58+o17*o3*o35*o58+o12*o22*o5*o58+o14*o3*o5*o58-o3*o51*o6-o17*o58*o6+o14*o22*o24*o65-o17*o24*o3*o65-o14*o3*o4*o65+o19*o5*o65+o17*o4*o5*o65-o22*o6*o65+o44*o9);
        x[3]=(-(o149*(-(o129*o167)+o133*(-(o4*o44)+o22*o58)))+o155*(-(o144*o167)+o133*(-(o44*o5)+o22*o65)))/(-(o149*(-(o127*o129)+o133*(o135-o3*o4)))+o155*(-(o127*o144)+o133*(o22*o35-o3*o5)));
        return;
    }

    template<class T> void linear_solve_5(std::vector<T> const &a,
                                          std::vector<T> const &b,
                                          std::vector<T> &x)
    {
        T o5=a[11];
        T o6=a[15];
        T o3=a[3];
        T o9=a[8];
        T o12=a[12];
        T o8=a[2];
        T o11=a[6];
        T o14=a[1];
        T o4=a[7];
        T o16=a[13];
        T o19=a[10];
        T o20=a[16];
        T o23=a[5];
        T o25=a[0];
        T o29=a[17];
        T o36=a[18];
        T o45=b[0];
        T o52=b[1];
        T o59=b[2];
        T o66=b[3];
        T o7=-(o3*o4*o5*o6);
        T o10=o5*o6*o8*o9;
        T o13=o11*o12*o3*o6;
        T o15=-(o12*o14*o6*o9);
        T o17=-(o11*o16*o6*o8);
        T o18=o14*o16*o4*o6;
        T o21=o19*o20*o3*o4;
        T o22=-(o19*o20*o8*o9);
        T o24=-(o12*o20*o23*o3);
        T o26=o12*o20*o25*o9;
        T o27=o16*o20*o23*o8;
        T o28=-(o16*o20*o25*o4);
        T o30=-(o11*o19*o29*o3);
        T o31=o14*o19*o29*o9;
        T o32=o23*o29*o3*o5;
        T o33=-(o25*o29*o5*o9);
        T o34=-(o14*o16*o23*o29);
        T o35=o11*o16*o25*o29;
        T o37=o11*o19*o36*o8;
        T o38=-(o14*o19*o36*o4);
        T o39=-(o23*o36*o5*o8);
        T o40=o25*o36*o4*o5;
        T o41=o12*o14*o23*o36;
        T o42=-(o11*o12*o25*o36);
        T o43=o10+o13+o15+o17+o18+o21+o22+o24+o26+o27+o28+o30+o31+o32+o33+o34+o35+o37+o38+o39+o40+o41+o42+o7;
        T o44=1/o43;
        T o75=a[4];
        T o77=a[9];
        T o81=a[14];
        T o96=a[19];
        T o111=-(o14*o23);
        T o112=o11*o25;
        T o113=o111+o112;
        T o119=-(o23*o8);
        T o120=o25*o4;
        T o121=o119+o120;
        T o107=-(o14*o19);
        T o108=o25*o5;
        T o109=o107+o108;
        T o104=-(o23*o75);
        T o105=o25*o77;
        T o106=o104+o105;
        T o122=-(o14*o6);
        T o123=o20*o25;
        T o124=o122+o123;
        T o155=a[20];
        T o132=-(o109*o121);
        T o133=-(o19*o8);
        T o134=o12*o25;
        T o135=o133+o134;
        T o136=o113*o135;
        T o137=o132+o136;
        T o146=-(o23*o3);
        T o147=o25*o9;
        T o148=o146+o147;
        T o156=-(o14*o155);
        T o157=a[21];
        T o158=o157*o25;
        T o159=o156+o158;
        T o149=-(o109*o148);
        T o150=-(o19*o3);
        T o151=o16*o25;
        T o152=o150+o151;
        T o153=o113*o152;
        T o154=o149+o153;
        T o125=-(o121*o124);
        T o126=-(o6*o8);
        T o127=o25*o29;
        T o128=o126+o127;
        T o129=o113*o128;
        T o130=o125+o129;
        T o110=-(o106*o109);
        T o114=-(o19*o75);
        T o115=o25*o81;
        T o116=o114+o115;
        T o117=o113*o116;
        T o118=o110+o117;
        T o160=-(o121*o159);
        T o161=-(o155*o8);
        T o162=a[22];
        T o163=o162*o25;
        T o164=o161+o163;
        T o165=o113*o164;
        T o166=o160+o165;
        T o167=-(o154*o166);
        T o168=-(o148*o159);
        T o169=-(o155*o3);
        T o170=a[23];
        T o171=o170*o25;
        T o172=o169+o171;
        T o173=o113*o172;
        T o174=o168+o173;
        T o175=o137*o174;
        T o176=o167+o175;
        T o200=-(o23*o45);
        T o201=o25*o52;
        T o202=o200+o201;
        T o178=-(o130*o154);
        T o179=-(o124*o148);
        T o180=-(o3*o6);
        T o181=o25*o36;
        T o182=o180+o181;
        T o183=o113*o182;
        T o184=o179+o183;
        T o185=o137*o184;
        T o186=o178+o185;
        T o203=-(o109*o202);
        T o204=-(o19*o45);
        T o205=o25*o59;
        T o206=o204+o205;
        T o207=o113*o206;
        T o208=o203+o207;
        T o131=-(o118*o130);
        T o138=-(o106*o124);
        T o139=-(o6*o75);
        T o140=o25*o96;
        T o141=o139+o140;
        T o142=o113*o141;
        T o143=o138+o142;
        T o144=o137*o143;
        T o145=o131+o144;
        T o177=-(o145*o176);
        T o187=-(o118*o166);
        T o188=-(o106*o159);
        T o189=-(o155*o75);
        T o190=a[24];
        T o191=o190*o25;
        T o192=o189+o191;
        T o193=o113*o192;
        T o194=o188+o193;
        T o195=o137*o194;
        T o196=o187+o195;
        T o197=o186*o196;
        T o198=o177+o197;
        T o199=1/o198;
        T o209=-(o130*o208);
        T o210=-(o124*o202);
        T o211=-(o45*o6);
        T o212=o25*o66;
        T o213=o211+o212;
        T o214=o113*o213;
        T o215=o210+o214;
        T o216=o137*o215;
        T o217=o209+o216;
        T o218=-(o176*o217);
        T o219=-(o166*o208);
        T o220=-(o159*o202);
        T o221=-(o155*o45);
        T o222=b[4];
        T o223=o222*o25;
        T o224=o221+o223;
        T o225=o113*o224;
        T o226=o220+o225;
        T o227=o137*o226;
        T o228=o219+o227;
        T o229=o186*o228;
        T o230=o218+o229;
        T o339=1/o186;
        x[0]=o44*(o11*o16*o29*o45-o11*o12*o36*o45-o16*o20*o4*o45+o36*o4*o45*o5-o14*o16*o29*o52-o12*o20*o3*o52+o12*o14*o36*o52+o29*o3*o5*o52-o11*o29*o3*o59+o20*o3*o4*o59-o14*o36*o4*o59+o11*o12*o3*o66+o14*o16*o4*o66-o3*o4*o5*o66+o16*o20*o52*o8-o36*o5*o52*o8+o11*o36*o59*o8-o11*o16*o66*o8+o12*o20*o45*o9-o29*o45*o5*o9+o14*o29*o59*o9-o12*o14*o66*o9-o20*o59*o8*o9+o5*o66*o8*o9)-o199*o230*o44*(o11*o16*o29*o75-o11*o12*o36*o75-o16*o20*o4*o75+o36*o4*o5*o75-o14*o16*o29*o77-o12*o20*o3*o77+o12*o14*o36*o77+o29*o3*o5*o77+o16*o20*o77*o8-o36*o5*o77*o8-o11*o29*o3*o81+o20*o3*o4*o81-o14*o36*o4*o81+o11*o36*o8*o81+o12*o20*o75*o9-o29*o5*o75*o9+o14*o29*o81*o9-o20*o8*o81*o9+o11*o12*o3*o96+o14*o16*o4*o96-o3*o4*o5*o96-o11*o16*o8*o96-o12*o14*o9*o96+o5*o8*o9*o96);
        x[1]=o44*(-(o16*o23*o29*o45)+o12*o23*o36*o45-o19*o36*o4*o45+o16*o25*o29*o52-o19*o29*o3*o52-o12*o25*o36*o52+o23*o29*o3*o59+o25*o36*o4*o59+o16*o4*o45*o6+o12*o3*o52*o6-o3*o4*o59*o6-o12*o23*o3*o66-o16*o25*o4*o66+o19*o3*o4*o66+o19*o36*o52*o8-o23*o36*o59*o8-o16*o52*o6*o8+o16*o23*o66*o8+o19*o29*o45*o9-o25*o29*o59*o9-o12*o45*o6*o9+o12*o25*o66*o9+o59*o6*o8*o9-o19*o66*o8*o9)-o199*o230*o44*(-(o16*o23*o29*o75)+o12*o23*o36*o75-o19*o36*o4*o75+o16*o4*o6*o75+o16*o25*o29*o77-o19*o29*o3*o77-o12*o25*o36*o77+o12*o3*o6*o77+o19*o36*o77*o8-o16*o6*o77*o8+o23*o29*o3*o81+o25*o36*o4*o81-o3*o4*o6*o81-o23*o36*o8*o81+o19*o29*o75*o9-o12*o6*o75*o9-o25*o29*o81*o9+o6*o8*o81*o9-o12*o23*o3*o96-o16*o25*o4*o96+o19*o3*o4*o96+o16*o23*o8*o96+o12*o25*o9*o96-o19*o8*o9*o96);
        x[2]=o44*(o16*o20*o23*o45+o11*o19*o36*o45-o23*o36*o45*o5-o16*o20*o25*o52+o19*o20*o3*o52-o14*o19*o36*o52+o25*o36*o5*o52-o20*o23*o3*o59+o14*o23*o36*o59-o11*o25*o36*o59-o11*o16*o45*o6+o14*o16*o52*o6-o3*o5*o52*o6+o11*o3*o59*o6-o14*o16*o23*o66+o11*o16*o25*o66-o11*o19*o3*o66+o23*o3*o5*o66-o19*o20*o45*o9+o20*o25*o59*o9+o45*o5*o6*o9-o14*o59*o6*o9+o14*o19*o66*o9-o25*o5*o66*o9)-o199*o230*o44*(o16*o20*o23*o75+o11*o19*o36*o75-o23*o36*o5*o75-o11*o16*o6*o75-o16*o20*o25*o77+o19*o20*o3*o77-o14*o19*o36*o77+o25*o36*o5*o77+o14*o16*o6*o77-o3*o5*o6*o77-o20*o23*o3*o81+o14*o23*o36*o81-o11*o25*o36*o81+o11*o3*o6*o81-o19*o20*o75*o9+o5*o6*o75*o9+o20*o25*o81*o9-o14*o6*o81*o9-o14*o16*o23*o96+o11*o16*o25*o96-o11*o19*o3*o96+o23*o3*o5*o96+o14*o19*o9*o96-o25*o5*o9*o96);
        x[3]=o217*o339-o145*o199*o230*o339;
        x[4]=o199*o230;
        return;
    }

    template<class T> void symmetric_linear_solve_5(std::vector<T> const &a,
                                                    std::vector<T> const &b,
                                                    std::vector<T> &x)
    {
        T o5=a[8];
        T o6=a[9];
        T o3=a[3];
        T o9=a[7];
        T o12=o6*o6;
        T o4=a[6];
        T o8=a[2];
        T o16=a[10];
        T o11=a[5];
        T o14=a[1];
        T o20=a[4];
        T o24=a[0];
        T o26=o16*o16;
        T o18=o9*o9;
        T o29=a[11];
        T o36=a[12];
        T o45=b[0];
        T o52=b[1];
        T o59=b[2];
        T o66=b[3];
        T o7=-(o3*o4*o5*o6);
        T o10=o5*o6*o8*o9;
        T o13=o11*o12*o3;
        T o15=-(o12*o14*o9);
        T o17=o16*o3*o4*o9;
        T o19=-(o16*o18*o8);
        T o21=-(o16*o20*o3*o6);
        T o22=-(o11*o16*o6*o8);
        T o23=o14*o16*o4*o6;
        T o25=o16*o24*o6*o9;
        T o27=o20*o26*o8;
        T o28=-(o24*o26*o4);
        T o30=-(o11*o29*o3*o9);
        T o31=o14*o18*o29;
        T o32=o20*o29*o3*o5;
        T o33=-(o24*o29*o5*o9);
        T o34=-(o14*o16*o20*o29);
        T o35=o11*o16*o24*o29;
        T o37=o11*o36*o8*o9;
        T o38=-(o14*o36*o4*o9);
        T o39=-(o20*o36*o5*o8);
        T o40=o24*o36*o4*o5;
        T o41=o14*o20*o36*o6;
        T o42=-(o11*o24*o36*o6);
        T o43=o10+o13+o15+o17+o19+o21+o22+o23+o25+o27+o28+o30+o31+o32+o33+o34+o35+o37+o38+o39+o40+o41+o42+o7;
        T o44=1/o43;
        T o86=o29*o29;
        T o80=o5*o5;
        T o95=a[13];
        T o105=o24*o5;
        T o110=-(o14*o20);
        T o111=o11*o24;
        T o112=o110+o111;
        T o114=o24*o29;
        T o118=-(o20*o8);
        T o119=o24*o4;
        T o120=o118+o119;
        T o107=-(o14*o9);
        T o108=o105+o107;
        T o103=o20*o20;
        T o104=-o103;
        T o106=o104+o105;
        T o121=-(o14*o6);
        T o122=o16*o24;
        T o123=o121+o122;
        T o130=-(o108*o120);
        T o131=-(o8*o9);
        T o132=o24*o6;
        T o133=o131+o132;
        T o134=o112*o133;
        T o135=o130+o134;
        T o144=-(o20*o3);
        T o145=o24*o9;
        T o146=o144+o145;
        T o152=-(o14*o16);
        T o153=o114+o152;
        T o138=o24*o95;
        T o147=-(o108*o146);
        T o148=-(o3*o9);
        T o149=o122+o148;
        T o150=o112*o149;
        T o151=o147+o150;
        T o124=-(o120*o123);
        T o125=-(o6*o8);
        T o126=o114+o125;
        T o127=o112*o126;
        T o128=o124+o127;
        T o156=o24*o36;
        T o109=-(o106*o108);
        T o113=-(o20*o9);
        T o115=o113+o114;
        T o116=o112*o115;
        T o117=o109+o116;
        T o154=-(o120*o153);
        T o155=-(o16*o8);
        T o157=o155+o156;
        T o158=o112*o157;
        T o159=o154+o158;
        T o160=-(o151*o159);
        T o161=-(o146*o153);
        T o162=-(o16*o3);
        T o163=o138+o162;
        T o164=o112*o163;
        T o165=o161+o164;
        T o166=o135*o165;
        T o167=o160+o166;
        T o190=-(o20*o45);
        T o191=o24*o52;
        T o192=o190+o191;
        T o169=-(o128*o151);
        T o170=-(o123*o146);
        T o171=-(o3*o6);
        T o172=o156+o171;
        T o173=o112*o172;
        T o174=o170+o173;
        T o175=o135*o174;
        T o176=o169+o175;
        T o193=-(o108*o192);
        T o194=-(o45*o9);
        T o195=o24*o59;
        T o196=o194+o195;
        T o197=o112*o196;
        T o198=o193+o197;
        T o129=-(o117*o128);
        T o136=-(o106*o123);
        T o137=-(o20*o6);
        T o139=o137+o138;
        T o140=o112*o139;
        T o141=o136+o140;
        T o142=o135*o141;
        T o143=o129+o142;
        T o168=-(o143*o167);
        T o177=-(o117*o159);
        T o178=-(o106*o153);
        T o179=-(o16*o20);
        T o180=a[14];
        T o181=o180*o24;
        T o182=o179+o181;
        T o183=o112*o182;
        T o184=o178+o183;
        T o185=o135*o184;
        T o186=o177+o185;
        T o187=o176*o186;
        T o188=o168+o187;
        T o189=1/o188;
        T o199=-(o128*o198);
        T o200=-(o123*o192);
        T o201=-(o45*o6);
        T o202=o24*o66;
        T o203=o201+o202;
        T o204=o112*o203;
        T o205=o200+o204;
        T o206=o135*o205;
        T o207=o199+o206;
        T o208=-(o167*o207);
        T o209=-(o159*o198);
        T o210=-(o153*o192);
        T o211=-(o16*o45);
        T o212=b[4];
        T o213=o212*o24;
        T o214=o211+o213;
        T o215=o112*o214;
        T o216=o210+o215;
        T o217=o135*o216;
        T o218=o209+o217;
        T o219=o176*o218;
        T o220=o208+o219;
        T o331=pow(o16,3.);
        T o352=pow(o29,3.);
        T o377=o36*o36;
        T o416=o95*o95;
        x[0]=o44*(o11*o16*o29*o45-o26*o4*o45+o36*o4*o45*o5-o14*o16*o29*o52+o29*o3*o5*o52-o11*o29*o3*o59+o16*o3*o4*o59-o14*o36*o4*o59-o11*o36*o45*o6-o16*o3*o52*o6+o14*o36*o52*o6+o14*o16*o4*o66-o3*o4*o5*o66+o11*o3*o6*o66+o26*o52*o8-o36*o5*o52*o8+o11*o36*o59*o8-o11*o16*o66*o8-o29*o45*o5*o9+o14*o29*o59*o9+o16*o45*o6*o9-o14*o6*o66*o9-o16*o59*o8*o9+o5*o66*o8*o9)-o189*o220*o44*(o11*o16*o20*o29-o20*o26*o4+o16*o29*o3*o4-o14*o29*o36*o4-o14*o16*o29*o5+o20*o36*o4*o5-o11*o20*o36*o6-o16*o3*o5*o6+o14*o36*o5*o6+o11*o29*o36*o8+o26*o5*o8+o29*o3*o80-o36*o8*o80-o11*o3*o86-o20*o29*o5*o9+o16*o20*o6*o9-o16*o29*o8*o9+o14*o86*o9+o14*o16*o4*o95-o3*o4*o5*o95+o11*o3*o6*o95-o11*o16*o8*o95-o14*o6*o9*o95+o5*o8*o9*o95);
        x[1]=o44*(o18*o29*o45-o16*o20*o29*o45+o16*o24*o29*o52+o12*o3*o52+o20*o29*o3*o59+o24*o36*o4*o59+o20*o36*o45*o6+o16*o4*o45*o6-o24*o36*o52*o6-o3*o4*o59*o6-o16*o24*o4*o66-o20*o3*o6*o66-o20*o36*o59*o8-o16*o52*o6*o8-o18*o66*o8+o16*o20*o66*o8-o12*o45*o9-o36*o4*o45*o9-o29*o3*o52*o9-o24*o29*o59*o9+o3*o4*o66*o9+o24*o6*o66*o9+o36*o52*o8*o9+o59*o6*o8*o9)-o189*o220*o44*(-(o103*o16*o29)+o18*o20*o29+o24*o29*o36*o4+o16*o24*o29*o5+o12*o3*o5+o103*o36*o6+o16*o20*o4*o6-o29*o3*o4*o6-o24*o36*o5*o6-o20*o29*o36*o8-o16*o5*o6*o8+o20*o3*o86-o12*o20*o9-o20*o36*o4*o9-o29*o3*o5*o9+o36*o5*o8*o9+o29*o6*o8*o9-o24*o86*o9-o16*o24*o4*o95-o20*o3*o6*o95-o18*o8*o95+o16*o20*o8*o95+o3*o4*o9*o95+o24*o6*o9*o95);
        x[2]=o44*(-(o16*o18*o45)+o20*o26*o45-o20*o36*o45*o5-o24*o26*o52+o24*o36*o5*o52-o16*o20*o3*o59+o14*o20*o36*o59-o11*o24*o36*o59-o11*o16*o45*o6+o14*o16*o52*o6-o3*o5*o52*o6+o11*o3*o59*o6+o14*o18*o66-o14*o16*o20*o66+o11*o16*o24*o66+o20*o3*o5*o66+o11*o36*o45*o9+o16*o3*o52*o9-o14*o36*o52*o9+o16*o24*o59*o9+o45*o5*o6*o9-o14*o59*o6*o9-o11*o3*o66*o9-o24*o5*o66*o9)-o189*o220*o44*(-(o16*o18*o20)+o103*o26-o16*o20*o29*o3+o14*o20*o29*o36-o11*o24*o29*o36-o24*o26*o5-o103*o36*o5-o11*o16*o20*o6+o11*o29*o3*o6+o14*o16*o5*o6+o24*o36*o80-o3*o6*o80+o16*o24*o29*o9+o11*o20*o36*o9+o16*o3*o5*o9-o14*o36*o5*o9-o14*o29*o6*o9+o20*o5*o6*o9+o14*o18*o95-o14*o16*o20*o95+o11*o16*o24*o95+o20*o3*o5*o95-o11*o3*o9*o95-o24*o5*o9*o95);
        x[3]=(-(o11*o12*o20*o212)+o16*o212*o24*o29*o4+o11*o12*o180*o45-o20*o352*o45+o16*o20*o29*o36*o45-o26*o29*o4*o45+o12*o14*o212*o5-o103*o212*o29*o5-o12*o29*o45*o5+o180*o20*o29*o45*o5-o12*o14*o180*o52+o12*o20*o29*o52+o24*o352*o52-o16*o24*o29*o36*o52+o16*o20*o29*o5*o52-o180*o24*o29*o5*o52-o11*o16*o20*o29*o59-o14*o180*o20*o29*o59+o11*o180*o24*o29*o59-o103*o16*o36*o59-o16*o180*o24*o4*o59+o20*o26*o4*o59+o14*o16*o29*o5*o59+o16*o24*o36*o5*o59+o103*o16*o212*o6-o14*o212*o29*o4*o6-o16*o180*o20*o45*o6-o11*o29*o36*o45*o6-o16*o212*o24*o5*o6+o20*o212*o4*o5*o6+o26*o45*o5*o6-o180*o4*o45*o5*o6+o16*o180*o24*o52*o6-o20*o26*o52*o6+o14*o29*o36*o52*o6-o20*o36*o5*o52*o6+o11*o20*o36*o59*o6+o14*o180*o4*o59*o6-o20*o29*o4*o59*o6-o14*o36*o5*o59*o6-o14*o20*o29*o36*o66+o11*o24*o29*o36*o66+o14*o16*o29*o4*o66+o103*o36*o5*o66-o16*o20*o4*o5*o66+o180*o24*o4*o5*o66+o11*o16*o20*o6*o66+o14*o180*o20*o6*o66-o11*o180*o24*o6*o66-o103*o29*o6*o66-o14*o16*o5*o6*o66+o24*o29*o5*o6*o66-o16*o20*o212*o29*o8+o26*o29*o52*o8+o16*o180*o20*o59*o8-o26*o5*o59*o8+o11*o212*o29*o6*o8+o180*o5*o52*o6*o8-o11*o180*o59*o6*o8+o29*o5*o59*o6*o8-o11*o16*o29*o66*o8-o180*o20*o5*o66*o8+o212*o24*o29*o80-o16*o29*o45*o80+o36*o45*o6*o80-o24*o36*o66*o80-o212*o6*o8*o80+o16*o66*o8*o80+o14*o20*o212*o86-o11*o212*o24*o86+o11*o16*o45*o86-o14*o16*o52*o86+o103*o59*o86-o24*o5*o59*o86+o4*o45*o6*o86-o24*o4*o66*o86-o52*o6*o8*o86+o20*o66*o8*o86+o11*o20*o212*o29*o9-o16*o20*o212*o4*o9-o11*o180*o29*o45*o9+o16*o180*o4*o45*o9-o14*o212*o29*o5*o9-o16*o36*o45*o5*o9+o14*o180*o29*o52*o9+o16*o20*o36*o52*o9-o11*o20*o36*o66*o9-o14*o180*o4*o66*o9+o20*o29*o4*o66*o9+o14*o36*o5*o66*o9+o16*o212*o5*o8*o9-o16*o180*o52*o8*o9+o11*o180*o66*o8*o9-o29*o5*o66*o8*o9+o45*o5*o86*o9-o20*o52*o86*o9-o212*o24*o4*o5*o95-o20*o36*o45*o5*o95+o16*o4*o45*o5*o95+o24*o36*o5*o52*o95+o14*o20*o36*o59*o95-o11*o24*o36*o59*o95-o14*o16*o4*o59*o95+o24*o29*o4*o59*o95-o14*o20*o212*o6*o95+o11*o212*o24*o6*o95-o11*o16*o45*o6*o95+o20*o29*o45*o6*o95+o14*o16*o52*o6*o95-o24*o29*o52*o6*o95+o20*o212*o5*o8*o95-o16*o5*o52*o8*o95+o11*o16*o59*o8*o95-o20*o29*o59*o8*o95+o14*o212*o4*o9*o95+o11*o36*o45*o9*o95-o29*o4*o45*o9*o95-o14*o36*o52*o9*o95-o11*o212*o8*o9*o95+o29*o52*o8*o9*o95)/(o14*o18*o180*o29-o14*o16*o180*o20*o29+o11*o16*o180*o24*o29-o11*o20*o26*o29+o11*o12*o180*o3-o20*o3*o352+o16*o18*o20*o36-o103*o26*o36+o16*o20*o29*o3*o36-o14*o20*o29*o377+o11*o24*o29*o377-o180*o24*o26*o4-o26*o29*o3*o4+o20*o331*o4+o14*o16*o29*o36*o4+o14*o26*o29*o5-o12*o29*o3*o5+o180*o20*o29*o3*o5+o24*o26*o36*o5+o103*o377*o5-o16*o20*o36*o4*o5+o180*o24*o36*o4*o5-o24*o4*o416*o5-o16*o180*o20*o3*o6+2.*o11*o16*o20*o36*o6+o14*o180*o20*o36*o6-o11*o180*o24*o36*o6-o103*o29*o36*o6-o11*o29*o3*o36*o6+o14*o16*o180*o4*o6-o16*o20*o29*o4*o6-o14*o20*o416*o6+o11*o24*o416*o6+o26*o3*o5*o6-2.*o14*o16*o36*o5*o6+o24*o29*o36*o5*o6-o180*o3*o4*o5*o6-o16*o18*o180*o8+o180*o20*o26*o8-o11*o16*o29*o36*o8-o331*o5*o8-o180*o20*o36*o5*o8+o20*o416*o5*o8-o11*o16*o180*o6*o8+o16*o29*o5*o6*o8-o16*o29*o3*o80-o24*o377*o80+o3*o36*o6*o80+o16*o36*o8*o80+o103*o16*o86-o18*o20*o86+o11*o16*o3*o86-o24*o36*o4*o86-o16*o24*o5*o86+o3*o4*o6*o86+o20*o36*o8*o86-o12*o14*o180*o9+o12*o20*o29*o9-o11*o180*o29*o3*o9+o24*o352*o9-o16*o24*o29*o36*o9-o11*o20*o377*o9+o16*o180*o3*o4*o9-o14*o180*o36*o4*o9+o20*o29*o36*o4*o9+o14*o4*o416*o9+o16*o20*o29*o5*o9-o180*o24*o29*o5*o9-o16*o3*o36*o5*o9+o14*o377*o5*o9+o16*o180*o24*o6*o9-o20*o26*o6*o9+o14*o29*o36*o6*o9-o20*o36*o5*o6*o9+o26*o29*o8*o9+o11*o180*o36*o8*o9-o11*o416*o8*o9-o29*o36*o5*o8*o9+o180*o5*o6*o8*o9-o14*o16*o86*o9+o3*o5*o86*o9-o6*o8*o86*o9-o11*o12*o20*o95-o14*o18*o36*o95+o14*o16*o20*o36*o95-o11*o16*o24*o36*o95-o14*o26*o4*o95+2.*o16*o24*o29*o4*o95+o12*o14*o5*o95-o103*o29*o5*o95-o20*o3*o36*o5*o95+o16*o3*o4*o5*o95+o103*o16*o6*o95-o11*o16*o3*o6*o95+o20*o29*o3*o6*o95-o14*o29*o4*o6*o95-o16*o24*o5*o6*o95+o20*o4*o5*o6*o95+o11*o26*o8*o95+o18*o29*o8*o95-2.*o16*o20*o29*o8*o95+o11*o29*o6*o8*o95+o24*o29*o80*o95-o6*o8*o80*o95+o14*o20*o86*o95-o11*o24*o86*o95+o11*o20*o29*o9*o95+o11*o3*o36*o9*o95-o16*o20*o4*o9*o95-o29*o3*o4*o9*o95-o14*o29*o5*o9*o95+o24*o36*o5*o9*o95+o14*o16*o6*o9*o95-o24*o29*o6*o9*o95);
        x[4]=o189*o220;
        return;
    }

    template<class T> void linear_solve(std::vector<T> const &a,
                                        std::vector<T> const &b,
                                        std::vector<T> &x)
    {
        int size = x.size();
        
        switch(size)
        {
        case 1:
            return linear_solve_1(a, b, x);
        case 2:
            return linear_solve_2(a, b, x);
        case 3:
            return linear_solve_3(a, b, x);
        case 4:
            return linear_solve_4(a, b, x);
        case 5:
            return linear_solve_5(a, b, x);
        default:
            AssertMsg(false, "linear solve of that size not implemented");
            return;
        }
    }

    template<class T> void symmetric_linear_solve(std::vector<T> const &a,
                                                  std::vector<T> const &b,
                                                  std::vector<T> &x)
    {
        int size = x.size();

        switch(size)
        {
        case 1:
            return symmetric_linear_solve_1(a, b, x);
        case 2:
            return symmetric_linear_solve_2(a, b, x);
        case 3:
            return symmetric_linear_solve_3(a, b, x);
        case 4:
            return symmetric_linear_solve_4(a, b, x);
        case 5:
            return symmetric_linear_solve_5(a, b, x);
        default:
            AssertMsg(false, "symmetric linear solve of that size not implemented");
            return;
        }
    }
}

#endif
