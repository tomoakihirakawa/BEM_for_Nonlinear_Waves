#ifndef lib_multipole_expansion_H
#define lib_multipole_expansion_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <ankerl/unordered_dense.h>
#include <numbers>
#include <sstream>

#include "basic_arithmetic_array_operations.hpp"
#include "basic_linear_systems.hpp"
#include "lib_multipole_expansion_constants.hpp"
#include "lib_quadmath.hpp"
using namespace std::literals::complex_literals;
using cmplx = std::complex<double>;

inline bool bemNearUseRunsEnabled() {
  static const bool enabled = []() -> bool {
    if (const char* env = std::getenv("BEM_NEAR_USE_RUNS")) {
      // Disable only when explicitly set to "0".
      return !(env[0] == '0' && env[1] == '\0');
    }
    return true;
  }();
  return enabled;
}

inline const std::array<std::array<double, 31>, 31> sqrt_nm_nm = {
    {{1., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.7071067811865475, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.4082482904638631, 0.20412414523193154, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.2886751345948129, 0.09128709291752768, 0.03726779962499649, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.22360679774997896, 0.05270462766947299, 0.014085904245475275, 0.004980119205559973, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.18257418583505536, 0.03450327796711771, 0.007042952122737638, 0.0016600397351866577, 0.00052495065695726, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.1543033499620919, 0.024397501823713332, 0.0040662503039522215, 0.0007423923386456234, 0.00015827857841616382, 0.00004569108992776174, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.1336306209562122, 0.018184824186332698, 0.0025717224993681985, 0.0003877017543326962, 0.00006461695905544937, 0.000012672428274337013, 3.3868489186454308e-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.1178511301977579, 0.014085904245475275, 0.0017338549553676645, 0.00022383971222927231, 0.00003104098307414404, 4.7897276744570195e-6, 8.744806305356235e-7, 2.1862015763390587e-7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.10540925533894598, 0.011236664374387369, 0.0012260205965343742, 0.00013881949648441295, 0.000016592103373156565, 2.1420313347595755e-6, 3.0917559193401366e-7, 5.3023176577281006e-8, 1.2497682572615703e-8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.09534625892455924, 0.009174698042719672, 0.0008996531605990218, 0.00009087869293935404, 9.579455348914039e-6, 1.0710156673797877e-6, 1.2987972715583107e-7, 1.7674392192427002e-8, 2.867165019077962e-9, 6.411175885367803e-10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.08703882797784893, 0.00763381020690574, 0.0006800738654737899, 0.00006208196614828809, 5.866194404653597e-6, 5.808397976391274e-7, 6.12258905403645e-8, 7.023091305098204e-9, 9.066771887846482e-10, 1.3990332756368328e-10, 2.982748965711089e-11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.08006407690254357, 0.0064517471761563645, 0.0005267829510341783, 0.000043898579252848194, 3.764272115814879e-6, 3.3534801352299794e-7, 3.14082191407746e-8, 3.14082191407746e-9, 3.4269176584847724e-10, 4.218244040462945e-11, 6.2194615286459366e-12, 1.2695422683377336e-12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.07412493166611012, 0.005524946201098299, 0.0004164584894532388, 0.000031940908071889604, 2.509514743876586e-6, 2.0354852405013068e-7, 1.720299011446777e-8, 1.532564167533272e-9, 1.4612425993612895e-10, 1.5234507220052806e-11, 1.7954038938891261e-12, 2.539084536675467e-13, 4.9795544612191934e-14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.06900655593423542, 0.0047847437646934545, 0.00033499907004387234, 0.0000238073472372723, 1.7271664998763255e-6, 1.2873539007279484e-7, 9.93215097345444e-9, 8.003555336576621e-10, 6.813078749549884e-11, 6.2194615286459366e-12, 6.219461528645936e-13, 7.042153453631633e-14, 9.583157028764344e-15, 1.8110464480707933e-15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.06454972243679027, 0.004184137043778615, 0.00027352559530478956, 0.000018114675073333608, 1.221291144300784e-6, 8.427709566117844e-8, 5.989312389111183e-9, 4.4153796733665665e-10, 3.406539374774942e-11, 2.7814277522994136e-12, 2.4394735152773264e-13, 2.3473844845438778e-14, 2.5612064489493374e-15, 3.363028826288546e-16, 6.140022498994532e-17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.06063390625908324, 0.0036900620230837303, 0.00022625221914921153, 0.00001403156697619594, 8.839056362752623e-7, 5.6819606303568e-8, 3.746573890260891e-9, 2.549220642992589e-10, 1.802571203400791e-11, 1.3361546727538149e-12, 1.049782255353683e-13, 8.872279396506633e-15, 8.237704614677945e-16, 8.683303091354011e-17, 1.1027805953831061e-17, 1.9494590928908317e-18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.057166195047502935, 0.0032787061473526402, 0.00018929618767677276, 0.00001103997785862415, 6.528075759249468e-7, 3.929439846428888e-8, 2.4184030470526603e-9, 1.5295323857955533e-10, 9.99886600236366e-12, 6.803366586737007e-13, 4.859547561955005e-14, 3.684013499396719e-15, 3.0079844263490085e-16, 2.7012497569312447e-17, 2.7569514884577653e-18, 3.39357269271957e-19, 5.819929153796601e-20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.05407380704358751, 0.0029325639294907786, 0.0001599844784189123, 8.806855686093134e-6, 4.907868795155064e-7, 2.7785335616744933e-8, 1.6041870997851784e-9, 9.485757175887514e-11, 5.772847978055655e-12, 3.6365524063136487e-13, 2.3875136219024277e-14, 1.647540922935589e-15, 1.2080356161406092e-16, 9.55036010402298e-18, 8.312521502205675e-19, 8.230622741349724e-20, 9.837475773099464e-21, 1.6395792955165772e-21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.051298917604257706, 0.002638531611922098, 0.00013643522161876318, 7.112177573145887e-6, 3.7484467091183e-7, 2.0036290454429938e-8, 1.0898305269346695e-9, 6.054614038525942e-11, 3.4499365303551642e-12, 2.0258724866669537e-13, 1.2329067328628607e-14, 7.828965582648662e-16, 5.230947661270771e-17, 3.717472628672112e-18, 2.8511713531899113e-19, 2.4096796001085158e-20, 2.3187152763056684e-21, 2.6954517578975646e-22, 4.37260014755555e-23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1., 0.048795003647426664, 0.0023866416420646975, 0.00011729703744891826, 5.807068671424468e-6, 2.903534335712234e-7, 1.4702617627274514e-8, 7.562210510333867e-10, 3.963675305041057e-11, 2.1247529892807687e-12, 1.1696380255210221e-13, 6.64309758380035e-15, 3.914482791324331e-16, 2.4091976159852608e-17, 1.561650865454643e-18, 1.07764147796743e-19, 8.032265333695052e-21, 6.602481433136971e-22, 6.183790431507657e-23, 7.001763889559219e-24, 1.1070760764863385e-24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1.,
      0.046524210519923545,
      0.0021692025004706648,
      0.00010158221421941786,
      4.7886315014996605e-6,
      2.277719728643592e-7,
      1.095868415392416e-8,
      5.347290332617261e-10,
      2.6538156330298823e-11,
      1.3438117822757118e-12,
      6.967339525147562e-14,
      3.713604444635547e-15,
      2.0442719657779699e-16,
      1.1686324992580568e-17,
      6.98391498455592e-19,
      4.39945291113154e-20,
      2.9527194609348998e-21,
      2.142127842145906e-22,
      1.7150748828864954e-23,
      1.5656420019457741e-24,
      1.728962355617497e-25,
      2.667846834301295e-26,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0},
     {1.,
      0.044455422447438706,
      0.001980201902173989,
      0.00008855732124870858,
      3.9843822502982624e-6,
      1.8073524913764726e-7,
      8.283986562225028e-9,
      3.845744215316325e-10,
      1.8129012089060738e-11,
      8.702204277741714e-13,
      4.2666066753341936e-14,
      2.144050525774124e-15,
      1.108662133605148e-16,
      5.9260483737718005e-18,
      3.2922490965398895e-19,
      1.9135809181461728e-20,
      1.1732917402031786e-21,
      7.670048049273388e-23,
      5.423542987667864e-24,
      4.2350755557433e-25,
      3.7729051754030186e-26,
      4.068428280973127e-27,
      6.133386406572036e-28,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0},
     {1.,
      0.04256282653793743,
      0.0018148850216015693,
      0.00007766990876952864,
      3.3423807018425985e-6,
      1.449105814256371e-7,
      6.342563843542926e-9,
      2.808534476165007e-10,
      1.2610692071105653e-11,
      5.755967094246757e-13,
      2.677918248384888e-14,
      1.2737558160678573e-15,
      6.215291969096202e-17,
      3.123301730909286e-18,
      1.623727252284406e-19,
      8.78011141314415e-21,
      4.970759255166987e-22,
      2.9705968359446517e-23,
      1.893983366497944e-24,
      1.3069726911875188e-25,
      9.965573343492674e-27,
      8.673918239448955e-28,
      9.143112624912157e-29,
      1.3480784107615843e-29,
      0,
      0,
      0,
      0,
      0,
      0,
      0},
     {1.,
      0.04082482904638631,
      0.0016694514082354448,
      0.0000684984209857492,
      2.8248272710847575e-6,
      1.1729456811558224e-7,
      4.912928827659812e-9,
      2.0798086168523098e-10,
      8.917105878934236e-12,
      3.8806722232841114e-13,
      1.7183905434211925e-14,
      7.762897181293485e-16,
      3.588400491449811e-17,
      1.7029795092108926e-18,
      8.32954505267255e-20,
      4.217829092369243e-21,
      2.222991118867921e-22,
      1.2274415085370757e-23,
      7.158584250066375e-25,
      4.4567398861618955e-26,
      3.0047334182847955e-27,
      2.2395960591833703e-28,
      1.906470771641401e-29,
      1.966374459242337e-30,
      2.83821705842792e-31,
      0,
      0,
      0,
      0,
      0,
      0},
     {1.,
      0.03922322702763681,
      0.0015408338784034144,
      0.0000607173606927742,
      2.403823329189731e-6,
      9.577061382776717e-8,
      3.8462404083659755e-9,
      1.5598564626392325e-10,
      6.400168590000365e-12,
      2.6621191827802367e-13,
      1.1249506769605595e-14,
      4.841016930192229e-16,
      2.127020725267102e-17,
      9.569918674447674e-19,
      4.423696472146341e-20,
      2.1089145461846216e-21,
      1.0415186570356942e-22,
      5.3569939276208894e-24,
      2.888297556246701e-25,
      1.645760272484188e-26,
      1.0015778060949317e-27,
      6.604208479255898e-29,
      4.816614068384924e-30,
      4.013845056987437e-31,
      4.054595797754172e-32,
      5.734064367124909e-33,
      0,
      0,
      0,
      0,
      0},
     {1.,
      0.037742567804819854,
      0.0014265349750363762,
      0.00005407266868359784,
      2.058511316097547e-6,
      7.8824472197141e-8,
      3.0407200297531575e-9,
      1.183598666326814e-10,
      4.656806667940934e-12,
      1.8553161858873915e-13,
      7.499671179737064e-15,
      3.082345867679868e-16,
      1.2910526133844702e-17,
      5.525195122815191e-19,
      2.4229583453505744e-20,
      1.0923539830333294e-21,
      5.082090666891961e-23,
      2.4508017457278914e-24,
      1.2315742160799112e-25,
      6.49096605058148e-27,
      3.6172739585553436e-28,
      2.1540552956581018e-29,
      1.3904367144822869e-30,
      9.93169081773062e-32,
      8.109191595508343e-33,
      8.029296133819685e-34,
      1.113463035245063e-34,
      0,
      0,
      0,
      0},
     {1.,
      0.036369648372665396,
      0.0013245044730538606,
      0.000048364065160539534,
      1.7731126566244665e-6,
      6.535779964392945e-8,
      2.4256554521012607e-9,
      9.077787379191073e-11,
      3.431081122865768e-12,
      1.3119066562861356e-13,
      5.083535034852086e-15,
      2.0000901082862465e-16,
      8.006768412092335e-18,
      3.268749516376755e-19,
      1.3643496875171323e-20,
      5.838877642269419e-22,
      2.5704225560302658e-23,
      1.1683738891046663e-24,
      5.507767332981392e-26,
      2.706919965402724e-27,
      1.3959873821192436e-28,
      7.615735533053188e-30,
      4.441587159991229e-31,
      2.809106370346196e-32,
      1.9667678521559944e-33,
      1.5746745256446799e-34,
      1.5294591045728086e-35,
      2.0813302159207996e-36,
      0,
      0,
      0},
     {1.,
      0.03509312031717982,
      0.0012330465600513202,
      0.00004343221265200662,
      1.5355606044085026e-6,
      5.456377014455781e-8,
      1.9511964368108084e-9,
      7.03162386800378e-11,
      2.5573768756813884e-12,
      9.401104366027673e-14,
      3.498728762031909e-15,
      1.3205100753166277e-16,
      5.063924975940288e-18,
      1.9771305335356593e-19,
      7.87707659356798e-21,
      3.210456743489421e-22,
      1.3423594001082063e-23,
      5.776595112493402e-25,
      2.5680097603370388e-26,
      1.184534573216848e-27,
      5.699094622593048e-29,
      2.8784774673280942e-30,
      1.5386109254700724e-31,
      8.795653226764121e-33,
      5.454832567602014e-34,
      3.7463943886574394e-35,
      2.943445419132117e-36,
      2.8064650906287853e-37,
      3.75029672752657e-38,
      0,
      0},
     {1.,
      0.033903175181040524,
      0.0011507487481140638,
      0.00003914926743091642,
      1.3365339529667335e-6,
      4.584273638574386e-8,
      1.5817246762906214e-9,
      5.4968732346059175e-11,
      1.926653071686668e-12,
      6.820277942358616e-14,
      2.442051239803031e-15,
      8.858250880144039e-17,
      3.2607672101605292e-18,
      1.2203114585475613e-19,
      4.652393819724967e-21,
      1.8109418382397407e-22,
      7.214964983015853e-24,
      2.950418589035446e-25,
      1.2423503413532843e-26,
      5.406635882463521e-28,
      2.4424691239684485e-29,
      1.1513909869312377e-30,
      5.700234783577265e-32,
      2.9877348446894465e-33,
      1.67543850471236e-34,
      1.0196394030485391e-35,
      6.874407452974273e-37,
      5.303720494991511e-38,
      4.9673893642040876e-39,
      6.52250178038964e-40,
      0},
     {1.,
      0.032791291789197645,
      0.001076426888425463,
      0.000035411844786990176,
      1.1687650369369893e-6,
      3.874418370480919e-8,
      1.2914727901603063e-9,
      4.333898916301203e-11,
      1.4659631761914363e-12,
      5.0047157644424506e-14,
      1.7267909916697387e-15,
      6.030215715220695e-17,
      2.134673223455801e-18,
      7.672924452278737e-20,
      2.8054990321811047e-21,
      1.0455477577944701e-22,
      3.9803322886358565e-24,
      1.5516950790121656e-25,
      6.2117517067664216e-27,
      2.5616832286005083e-28,
      1.0923053988275624e-29,
      4.836809603787499e-31,
      2.2358160302708322e-32,
      1.0858082505169e-33,
      5.584795015707868e-35,
      3.074328473976332e-36,
      1.8372626732939303e-37,
      1.2167569296028167e-38,
      9.224210478429689e-40,
      8.49157403659412e-41,
      1.0962574942272106e-41}}};

using namespace std::literals::complex_literals;
using cmplx = std::complex<double>;

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT Multipole_Expansion

# 多重極展開

この実装は，\cite{Greengard1997a}に基づいている．

## Green関数の多重極展開

次のGreen関数を考える．

$$
G({\bf x},{\bf a}) = \frac{1}{\|{\bf x}-{\bf a}\|},
\quad \nabla G({\bf x},{\bf a}) = -\frac{{\bf x}-{\bf a}}{\|{\bf x}-{\bf a}\|^3}
$$

グリーン関数は，球面調和関数を使って近似できる．
近似を$G_{\rm apx}({\bf x},{\bf a},{\bf c})$とする．

$$
G_{\rm apx}(n, {\bf x},{\bf a},{\bf c}) = \sum_{k=0}^n \sum_{m=-k}^k \left( \frac{r_{\rm near}}{r_{\rm far}} \right)^k \frac{1}{r_{\rm far}} Y(k, -m, a_{\rm near}, b_{\rm near}) Y(k, m, a_{\rm far}, b_{\rm far})=
{\bf Y}^\ast({\bf x},{\bf c})\cdot{\bf Y}({\bf a},{\bf c})
$$

$$
{\bf Y}^\ast({\bf x},{\bf c}) = r_{\rm near}^k Y(k, -m, a_{\rm near},b_{\rm near}), \quad {\bf Y}({\bf a},{\bf c}) = r_{\rm far}^{-k-1} Y(k, m, a_{\rm far}, b_{\rm far})
$$

ここで，$(r_{\rm near},a_{\rm near},b_{\rm near})$は，球面座標系に${\bf x}-{\bf c}$を変換したものであり，
$(r_{\rm far},a_{\rm far},b_{\rm far})$は，球面座標系に${\bf a}-{\bf c}$を変換したもの．$Y(k, m, a, b)$は球面調和関数：

$$
Y(k, m, a, b) = \sqrt{\frac{(k - |m|)!}{(k + |m|)!}} P_k^{|m|}(\cos(a)) e^{i mb}
$$

$P_k^m(x)$はルジャンドル陪関数：

$$
P_k^m(x) = \frac{(-1)^m}{2^k k!} (1-x^2)^{m/2} \frac{d^{k+m}}{dx^{k+m}}(x^2-1)^k
$$

### 球面座標系への変換

${\bf x}=(x,y,z)$から球面座標$(r,a,b)$への変換は次のように行う．

$$
r = \|{\bf x}\|, \quad a = \arctan \frac{\sqrt{x^2 + y^2}}{z}, \quad b = \arctan \frac{y}{x}
$$

$r_\parallel=\sqrt{x^2+y^2}$とする．$\frac{\partial}{\partial t}(\arctan(f(t))) = \frac{f'(t)}{1 + f(t)^2}$なので，
$(r,a,b)$の$(x,y,z)$に関する勾配は次のようになる．

$$
\nabla r = \frac{\bf x}{r},\quad
\nabla a = \frac{1}{r^2r_\parallel} \left(xz,yz,-r_\parallel^2\right),\quad
\nabla b = \frac{1}{r_\parallel^2} \left(-y,x,0\right)
$$

*/

inline constexpr double factorial(const int n) {
  if (n < 0)
    throw std::runtime_error("factorial of negative number");
  else if (n < factorial_list_as_double.size())
    return factorial_list_as_double[n];
  else
    return n * factorial(n - 1);
}

#if defined(_LIBCPP_VERSION)
namespace bem_sf {
inline double assoc_legendre_nonneg_m(const int l, const int m, const double x) {
  if (m < 0 || m > l)
    return 0.0;
  // Numerical Recipes recurrence with Condon–Shortley phase.
  double pmm = 1.0;
  if (m > 0) {
    const double somx2 = std::sqrt(std::max(0.0, (1.0 - x) * (1.0 + x)));
    double fact = 1.0;
    for (int i = 1; i <= m; ++i) {
      pmm *= (-fact) * somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  double pmmp1 = x * (2.0 * m + 1.0) * pmm;
  if (l == m + 1)
    return pmmp1;
  double pll = 0.0;
  double pprev = pmm;
  double pcurr = pmmp1;
  for (int ll = m + 2; ll <= l; ++ll) {
    pll = ((2.0 * ll - 1.0) * x * pcurr - (ll + m - 1.0) * pprev) / (ll - m);
    pprev = pcurr;
    pcurr = pll;
  }
  return pcurr;
}

inline double assoc_legendre(const int l, const int m, const double x) {
  if (m >= 0)
    return assoc_legendre_nonneg_m(l, m, x);
  const int mp = -m;
  const double p = assoc_legendre_nonneg_m(l, mp, x);
  const double scale = (mp & 1 ? -1.0 : 1.0) * (factorial(l - mp) / factorial(l + mp));
  return scale * p;
}

inline double sph_legendre(const int l, const int m, const double theta) {
  const double x = std::cos(theta);
  const double p = assoc_legendre_nonneg_m(l, m, x);
  const double norm = std::sqrt((2.0 * l + 1.0) / (4.0 * std::numbers::pi) * (factorial(l - m) / factorial(l + m)));
  return norm * p;
}
} // namespace bem_sf
#endif

/* ------------------------- Basic Vector Operations ------------------------ */

/*

\cite{Liu_2009}p.70には，体球調和関数を使った高速多重極展開のアルゴリズムが書かれている．
以下は，その定義通りの体球調和関数を計算する関数である．

*/

inline constexpr cmplx complex_zero(0., 0.);
inline constexpr std::array<cmplx, 3> complex_zero3D{complex_zero, complex_zero, complex_zero};

inline cmplx Dot(const std::array<double, 3>& normal, const std::array<cmplx, 3>& abc) { return std::get<0>(abc) * std::get<0>(normal) + std::get<1>(abc) * std::get<1>(normal) + std::get<2>(abc) * std::get<2>(normal); };

inline cmplx Dot(const std::array<cmplx, 3>& abc, const std::array<double, 3>& normal) { return std::get<0>(abc) * std::get<0>(normal) + std::get<1>(abc) * std::get<1>(normal) + std::get<2>(abc) * std::get<2>(normal); };

inline std::array<cmplx, 3> operator*(std::array<cmplx, 3> a, const std::array<cmplx, 3>& A) {
  std::get<0>(a) *= std::get<0>(A);
  std::get<1>(a) *= std::get<1>(A);
  std::get<2>(a) *= std::get<2>(A);
  return a;
};

inline std::array<cmplx, 3> operator*(cmplx a, std::array<cmplx, 3> A) {
  std::get<0>(A) *= a;
  std::get<1>(A) *= a;
  std::get<2>(A) *= a;
  return A;
};

inline std::array<cmplx, 3> operator*(cmplx a, std::array<double, 3> A) {
  std::array<cmplx, 3> B = {a, a, a};
  std::get<0>(B) *= std::get<0>(A);
  std::get<1>(B) *= std::get<1>(A);
  std::get<2>(B) *= std::get<2>(A);
  return B;
};

inline std::array<cmplx, 3> operator+(std::array<cmplx, 3> a, const std::array<double, 3>& A) {
  std::get<0>(a) += std::get<0>(A);
  std::get<1>(a) += std::get<1>(A);
  std::get<2>(a) += std::get<2>(A);
  return a;
};

inline std::array<cmplx, 3> operator+(const std::array<double, 3>& A, const std::array<cmplx, 3>& a) { return a + A; };

// 2つの2次元ベクトルの要素ごとの積を計算する operator*
template <typename T>
inline std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& lhs, const std::vector<std::vector<T>>& rhs) {
  // --- サイズのチェック ---
  if (lhs.empty() || rhs.empty() || lhs.size() != rhs.size() || lhs[0].size() != rhs[0].size()) {
    throw std::invalid_argument("Matrices must have the same dimensions for element-wise multiplication.");
  }

  // 結果を格納する新しい2次元ベクトルを作成
  std::vector<std::vector<T>> result = lhs; // サイズをlhsと同じに初期化

  // --- 要素ごとの積 ---
  for (size_t i = 0; i < result.size(); ++i) {
    for (size_t j = 0; j < result[i].size(); ++j) {
      result[i][j] = lhs[i][j] * rhs[i][j];
    }
  }

  return result;
}

// 2つの2次元ベクトルを要素ごとに加算する operator+=
template <typename T>
std::vector<std::vector<T>>& operator+=(std::vector<std::vector<T>>& lhs, const std::vector<std::vector<T>>& rhs) {
  // --- サイズのチェック ---
  if (lhs.empty() || rhs.empty() || lhs.size() != rhs.size() || lhs[0].size() != rhs[0].size()) {
    throw std::invalid_argument("Matrices must have the same dimensions for element-wise addition.");
  }

  // --- 要素ごとの加算 ---
  for (size_t i = 0; i < lhs.size(); ++i) {
    for (size_t j = 0; j < lhs[i].size(); ++j) {
      lhs[i][j] += rhs[i][j];
    }
  }

  return lhs; // 左辺への参照を返す
}
/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT Multipole_Expansion

## C++上での，Greengardの球面調和関数

`sph_harmonics_`

Greengardｎ(1997)の(3.15)と同じように，球面調和関数を定義する．
c++の`std::sph_legendre`を使って(3.15)を使う場合，係数を調整と，mの絶対値を考慮する必要がある．

c++での球面調和関数の定義は次のようになる[球面調和関数](https://cpprefjp.github.io/reference/cmath/sph_legendre.html)．
ただし，$\phi=0$の結果が返ってくるので，$e^{im\phi}$をかける必要がある．

$$
\begin{align*}
{\mathrm{std::sph\_legendre(n,m,\theta)}} &= (-1)^m \sqrt{\frac{(2n+1)(n-m)!}{4\pi(n+m)!}} {\rm{std::assoc_legendre}(n,m,cos(\theta))}\\
& = (-1)^m \sqrt{\frac{(2n+1)(n-m)!}{4\pi(n+m)!}} (1-x^2)^{m/2} \frac{d^m}{dx^m} P_n(x), \quad x = \cos(\theta)
\end{align*}
$$

Greengardｎ(1997)の(3.15)：

$$
\begin{align*}
Y(n, m, \theta, \phi) &= \sqrt{\frac{(n-|m|)!}{(n+|m|)!}} P_n^{|m|}(\cos(\theta)) e^{im \phi}\\
& = (-1)^{|m|}\sqrt{\frac{(n-|m|)!}{(n+|m|)!}} (1-x^2)^{|m|/2} \frac{d^{|m|}}{dx^{|m|}} P_n(x) e^{im \phi}, \quad x = \cos(\theta)
\end{align*}
$$

従って，$Y(n, m, \theta, \phi)$はc++の`std::sph_legendre`を使って

$$
Y(n, m, \theta, \phi) = \sqrt{\frac{4\pi}{2n+1}}{\mathrm{std::sph\_legendre(n,|m|,\theta)}} e^{im\phi}
$$

と計算できる．

*/

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT

`source4FMM`は，境界要素法における**節点上の値を表すわけではない**ことに注意する．
境界要素法において，面積分を離散化し，和に落とし込んだときのある１項分がソースと呼ばれるべきものであり，`source4FMM`はこのソースを表す．

しかし，BEMでは，ソースではなく，節点上の値に注目して扱うことが多い．
なので，`source4FMM`は，ソースとしての自身の値を，節点上の値から効率的に計算する関数`set_source_densities`を持っている．

ラプラス方程式を満たすポテンシャルは，グリーンの定理を使うと，
単層ポテンシャル(single layer potential) と２重層ポテンシャル(double layer potential)の和で表される．
(ただし，立体角にも影響される)

* 単層ポテンシャルは，極またはソース(monosource)の積分であり，
* ２重層ポテンシャルは，双極子または無限に近い２つの極(disource)の積分である．

単層ポテンシャルも２重層ポテンシャルも，
原点からの距離が近いところでは，影響が大きく，その影響度が，ラプラス方程式の基本解またはその法線微分となっている(カーネル関数)．

*/

struct DirectAccumulator {
  std::vector<double> phi_acc, phin_acc;
  std::vector<uint8_t> touched_flags;
  std::vector<int32_t> touched_list;
  std::size_t high_water_mark = 0;

  inline void reset() {
    if (high_water_mark > 0 && phi_acc.size() < high_water_mark) {
      phi_acc.assign(high_water_mark, 0.);
      phin_acc.assign(high_water_mark, 0.);
      touched_flags.assign(high_water_mark, 0);
    }
    touched_list.clear();
  }

  inline void add(int32_t idx, const double& phi_w, const double& phin_w) {
    const std::size_t i = static_cast<std::size_t>(idx);
    if (i >= phi_acc.size()) [[unlikely]] {
      const std::size_t new_size = i + 1;
      phi_acc.resize(new_size, 0.);
      phin_acc.resize(new_size, 0.);
      touched_flags.resize(new_size, 0);
    }
    phi_acc[i] += phi_w;
    phin_acc[i] += phin_w;
    if (!touched_flags[i]) {
      touched_flags[i] = 1;
      touched_list.push_back(idx);
    }
  }

  inline void update_high_water_mark() {
    if (phi_acc.size() > high_water_mark)
      high_water_mark = phi_acc.size();
  }
};

template <typename T>
struct source4FMM {
  // --- ツリー配置用（代表位置 = 要素重心） ---
  Tddd X{0.0, 0.0, 0.0}; // public: Buckets が p->X にアクセス

  // --- 内部P2Mソース（複数求積点を内包） ---
  struct P2MSource {
    Tddd X{0.0, 0.0, 0.0};      // 求積点の座標
    Tddd normal{0.0, 0.0, 1.0};  // その点での法線
    Tdd weighted_source_densities{0.0, 0.0};
    std::function<std::array<double, 2>()> get_weighted_source_densities;

    inline void updateDensity() {
      if (get_weighted_source_densities)
        weighted_source_densities = get_weighted_source_densities();
    }
  };
  std::vector<P2MSource> p2m_sources; // D1=1個, D3=3個, D6=6個, D7=7個

  // --- MAC拡張用: 重心から最も遠い内部P2Mソースまでの距離 ---
  double max_source_offset = 0.0;

  // --- 直接積分（要素レベル、1面1回） ---
  std::function<void(const T*, DirectAccumulator&)> fill_direct_entries;
  std::function<void(const T*, DirectAccumulator&)> fill_direct_entries_nonadj;

  // --- 面情報（既存） ---
  std::array<const T*, 3> face_vertices{nullptr, nullptr, nullptr};
  std::array<const T*, 3> face_midpoints{nullptr, nullptr, nullptr}; // 辺中点 (true_quadratic のみ)
  int32_t dof_indices[3] = {-1, -1, -1}; // Linear element DOF indices (for SIMD path)
  float near_region = 0.f;                // near/far threshold (for SIMD path)

  // --- 後方互換用（旧APIとの橋渡し） ---
  Tddd normal{0.0, 0.0, 1.0};
  Tdd weighted_source_densities{0.0, 0.0};
  std::function<std::array<double, 3>()> get_X;
  std::function<std::array<double, 3>()> get_normal;
  std::function<std::array<double, 2>()> get_weighted_source_densities;

  [[nodiscard]] inline bool isAdjacentTo(const T* target) const noexcept {
    return face_vertices[0] == target || face_vertices[1] == target || face_vertices[2] == target ||
           face_midpoints[0] == target || face_midpoints[1] == target || face_midpoints[2] == target;
  }

  // --- Ctors ---
  source4FMM() = default;

  source4FMM(std::function<std::array<double, 3>()> get_X_,
             std::function<std::array<double, 3>()> get_normal_,
             std::function<std::array<double, 2>()> get_weighted_source_densities_,
             std::function<void(const T*, DirectAccumulator&)> direct_handler_)
      : get_X(std::move(get_X_)),
        get_normal(std::move(get_normal_)),
        get_weighted_source_densities(std::move(get_weighted_source_densities_)),
        fill_direct_entries(std::move(direct_handler_)) { update(); }

  // --- Mutators ---
  inline void set_geometry(const Tddd& X_, const Tddd& n_) noexcept {
    X = X_;
    normal = n_;
  }

  inline void set_callbacks(std::function<std::array<double, 2>()> get_weighted_source_densities_,
                            std::function<void(const T*, DirectAccumulator&)> direct_handler_,
                            std::function<void(const T*, DirectAccumulator&)> direct_handler_nonadj_ = nullptr) {
    get_weighted_source_densities = std::move(get_weighted_source_densities_);
    fill_direct_entries = std::move(direct_handler_);
    fill_direct_entries_nonadj = std::move(direct_handler_nonadj_);
  }

  // --- Lifecycle / helpers ---
  inline void update() {
    updateDensity();
    updatePosition();
  }

  inline void updateDensity() {
    if (!p2m_sources.empty()) {
      for (auto& s : p2m_sources) s.updateDensity();
    } else if (get_weighted_source_densities) {
      weighted_source_densities = get_weighted_source_densities();
    }
  }

  inline void updatePosition() {
    if (get_X) this->X = get_X();
    if (get_normal) this->normal = get_normal();
  }

  [[nodiscard]] inline bool valid() const noexcept {
    return !p2m_sources.empty() || static_cast<bool>(get_weighted_source_densities);
  }

  // 可読性の高いアクセサ（後方互換用）
  inline double monopole() const noexcept { return weighted_source_densities[0]; }
  inline double dipole_strength() const noexcept { return weighted_source_densities[1]; }
  inline Tddd dipole_vec() const noexcept {
    return Tddd{normal[0] * weighted_source_densities[1], normal[1] * weighted_source_densities[1], normal[2] * weighted_source_densities[1]};
  }

};

/* -------------------------------------------------------------------------- */

struct SphericalCoordinates {
  std::array<double, 3> X;
  double r2D;
  double rho;
  double div_r2D;
  double x_r, y_r, z_r;
  double phi;
  double theta;
  double scale;
  SphericalCoordinates(const std::array<double, 3>& XIN, const double scale = 1.)
      : X(XIN / scale), r2D(std::hypot(std::get<0>(X), std::get<1>(X))), //! 方向非依存
        rho(std::hypot(std::get<0>(X), std::get<1>(X), std::get<2>(X))), //! 方向非依存
        div_r2D(1. / r2D),                                               //! 方向非依存
        x_r(std::get<0>(X) / rho),                                       //% 方向依存(-1倍)
        y_r(std::get<1>(X) / rho),                                       //% 方向依存(-1倍)
        z_r(std::get<2>(X) / rho),                                       //% 方向依存(-1倍)
        phi(std::atan2(std::get<1>(X), std::get<0>(X))),                 //! 方向非依存
        theta(std::atan2(r2D, std::get<2>(X))),                          //% 方向依存(-1倍) これは元の定義z/r=cos(theta)を書き換えたもの
        scale(scale)                                                     //! 方向非依存
  // theta(x_r * std::cos(phi) >= 0 ? std::atan2(r2D, std::get<2>(X)) : std::atan2(-r2D, std::get<2>(X)))  //! これは元の定義z/r=cos(theta)を書き換えたもの
  {}

  void initialize(const std::array<double, 3>& XIN, const double scale = 1.) {
    this->X = XIN / scale;
    this->r2D = std::hypot(std::get<0>(X), std::get<1>(X));                 //! 方向非依存
    this->rho = std::hypot(std::get<0>(X), std::get<1>(X), std::get<2>(X)); //! 方向非依存
    this->div_r2D = 1. / r2D;                                               //! 方向非依存
    this->x_r = std::get<0>(X) / rho;                                       //% 方向依存(-1倍)
    this->y_r = std::get<1>(X) / rho;                                       //% 方向依存(-1倍)
    this->z_r = std::get<2>(X) / rho;                                       //% 方向依存(-1倍)
    this->phi = std::atan2(std::get<1>(X), std::get<0>(X));                 //! 方向非依存
    this->theta = std::atan2(r2D, std::get<2>(X));                          //% 方向依存(-1倍) これは元の定義z/r=cos(theta)を書き換えたもの
    this->scale = scale;                                                    //! 方向非依存
  }

  int pre_max = -1;     // 0,1,2,...,pre_max (for sph_harmonics, rho, div_rhon1)
  int pre_max_ypq = -1; // independent bound for Ypq (i_absk_A_FMM)
  unsigned pre_mask = 0;
  std::vector<std::vector<cmplx>> pre_compute_sph_harmonics, pre_compute_sph_harmonics_rho_n, pre_compute_sph_harmonics_div_rhon1, pre_compute_div_rhon1_i_absk_A_FMM;

  enum PrecomputeMask : unsigned {
    PrecomputeBase = 1u << 0,
    PrecomputeRho = 1u << 1,
    PrecomputeDiv = 1u << 2,
    PrecomputeYpq = 1u << 3
  };

  void precompute_sph_mask(const int pre_max_IN, unsigned mask) {
    if (pre_max_IN < 0)
      return;
    // Keep the larger bound if already computed.
    const int target_max = (pre_max < 0) ? pre_max_IN : std::max(pre_max, pre_max_IN);
    mask |= PrecomputeBase;
    pre_max = target_max;
    pre_mask |= mask;
    if (pre_mask & PrecomputeBase)
      pre_compute_sph_harmonics.assign(pre_max + 1, std::vector<cmplx>(2 * pre_max + 1, 0.));
    if (pre_mask & PrecomputeRho)
      pre_compute_sph_harmonics_rho_n.assign(pre_max + 1, std::vector<cmplx>(2 * pre_max + 1, 0.));
    if (pre_mask & PrecomputeDiv)
      pre_compute_sph_harmonics_div_rhon1.assign(pre_max + 1, std::vector<cmplx>(2 * pre_max + 1, 0.));
    if (pre_mask & PrecomputeYpq)
      pre_compute_div_rhon1_i_absk_A_FMM.assign(pre_max + 1, std::vector<cmplx>(2 * pre_max + 1, 0.));

    if (pre_mask & PrecomputeYpq)
      pre_max_ypq = std::min(pre_max, N_i_absk_A_FMM);

    for (int n = 0; n <= pre_max; n++) {
      const double rho_n = std::pow(rho, n); // runtime value, cannot be constexpr
      const double inv_rho_n1 = 1.0 / (rho_n * rho);
      for (int m = -n; m <= n; m++) {
        auto tmp = sph_harmonics_polar_form_(n, m);
        const cmplx Y = std::polar(std::get<0>(tmp), std::get<1>(tmp));
        if (pre_mask & PrecomputeBase)
          this->pre_compute_sph_harmonics[n][m + pre_max] = Y;
        if (pre_mask & PrecomputeRho)
          this->pre_compute_sph_harmonics_rho_n[n][m + pre_max] = Y * rho_n;
        if (pre_mask & PrecomputeDiv)
          this->pre_compute_sph_harmonics_div_rhon1[n][m + pre_max] = Y * inv_rho_n1;
        if ((pre_mask & PrecomputeYpq) && n <= pre_max_ypq)
          this->pre_compute_div_rhon1_i_absk_A_FMM[n][m + pre_max] = Y * inv_rho_n1 / static_cast<std::complex<double>>(i_absk_A_FMM[n][m + N_i_absk_A_FMM]);
      }
    }
  }

  void precompute_sph_rho(const int pre_max_IN) { precompute_sph_mask(pre_max_IN, PrecomputeRho); }
  void precompute_sph_div_rhon1(const int pre_max_IN) { precompute_sph_mask(pre_max_IN, PrecomputeDiv); }
  void precompute_sph_Ypq(const int pre_max_IN) { precompute_sph_mask(pre_max_IN, PrecomputeYpq); }
  void precompute_sph(const int pre_max_IN) { precompute_sph_mask(pre_max_IN, PrecomputeRho | PrecomputeDiv | PrecomputeYpq); }

private:
  std::array<std::array<double, 3>, 3> Jacobian_inv() const {
    return {{/*grad r*/ {x_r, y_r, z_r},
             /*grad theta*/ {x_r * z_r * div_r2D, y_r * z_r * div_r2D, -r2D / rho / rho},
             /*grad phi*/ {-std::get<1>(X) * div_r2D * div_r2D, std::get<0>(X) * div_r2D * div_r2D, 0.}}};
  }

  inline int negOnePower(const int p) const { return 1 - ((p & 1) << 1); }

  inline cmplx sph_harmonics_(const int n, const int m) const {
    const double sqrt_MPI2 = std::sqrt(2.0 * std::numbers::pi);
    const double abs_m = std::abs(m);
    if (abs_m <= n) {
      double s = sqrt_MPI2 / std::sqrt(n + 0.5);
#if defined(_LIBCPP_VERSION)
      return std::polar(bem_sf::sph_legendre(n, static_cast<int>(abs_m), theta) * s, m * phi);
#else
      return std::polar(std::sph_legendre(n, abs_m, theta) * s, m * phi);
#endif
    } else
      return 0.0;
  }

  inline std::array<double, 2> sph_harmonics_polar_form_(const int n, const int m) const {
    const double sqrt_MPI2 = std::sqrt(2.0 * std::numbers::pi);
    const double abs_m = std::abs(m);
    if (abs_m <= n) {
      double s = sqrt_MPI2 / std::sqrt(n + 0.5);
#if defined(_LIBCPP_VERSION)
      return {bem_sf::sph_legendre(n, static_cast<int>(abs_m), theta) * s, m * phi};
#else
      return {std::sph_legendre(n, abs_m, theta) * s, m * phi};
#endif
    } else
      return {0.0, 0.0};
  }

  inline cmplx sph_harmonics(const int n, const int m) const {
    if ((pre_mask & PrecomputeBase) && n <= pre_max && std::abs(m) <= pre_max)
      return this->pre_compute_sph_harmonics[n][m + pre_max];
    else
      return this->sph_harmonics_(n, m);
  }

  inline cmplx sph_harmonics_rho(const int n, const int m) const {
    if (!(pre_mask & PrecomputeRho) || n > pre_max || m > pre_max || m < -pre_max)
      return this->sph_harmonics_(n, m) * std::pow(rho, n);
    else
      return this->pre_compute_sph_harmonics_rho_n[n][m + pre_max];
  }

  inline cmplx sph_harmonics_div_rhon1(const int n, const int m) const {
    if ((pre_mask & PrecomputeDiv) && n <= pre_max && std::abs(m) <= pre_max)
      return this->pre_compute_sph_harmonics_div_rhon1[n][m + pre_max];
    else
      return this->sph_harmonics_(n, m) / (std::pow(rho, n) * rho);
  }

  inline std::array<cmplx, 2> SolidHarmonicR_ForNear_Grad_SolidHarmonicR_ForNear_normal(const int n, const int m, const std::array<double, 3>& normal) const {
    cmplx Rnm = sph_harmonics_rho(n, m);
    auto [J0, J1, J2] = Jacobian_inv();
    /*rについての微分*/
    cmplx grad_Rnm1_dot_normal = (static_cast<double>(n) * (Rnm / rho)) * Dot(normal, J0);
    /*thetaについての微分*/
    const double cos_theta = std::get<2>(X) / rho;
    const int abs_m = std::abs(m);
#if defined(_LIBCPP_VERSION)
    grad_Rnm1_dot_normal += (std::get<2>(X) / r2D * abs_m - (bem_sf::assoc_legendre(n, abs_m + 1, cos_theta) / bem_sf::assoc_legendre(n, abs_m, cos_theta))) * Rnm * Dot(normal, J1);
#else
    grad_Rnm1_dot_normal += (std::get<2>(X) / r2D * abs_m - (std::assoc_legendre(n, abs_m + 1, cos_theta) / std::assoc_legendre(n, abs_m, cos_theta))) * Rnm * Dot(normal, J1);
#endif
    /*phiについての微分*/
    grad_Rnm1_dot_normal += (cmplx(0., m) * Rnm) * Dot(normal, J2);
    return {Rnm, grad_Rnm1_dot_normal};
  }

  //! Nishimura2002 eq. (15)
  inline cmplx Rnm_Nishimura2002(const int n, const int m) const {
    const double cos_theta = std::get<2>(X) / rho;
#if defined(_LIBCPP_VERSION)
    return std::polar(bem_sf::assoc_legendre(n, m, cos_theta) * (std::pow(rho, n) / factorial(n + m)), m * phi);
#else
    return std::polar(std::assoc_legendre(n, m, cos_theta) * (std::pow(rho, n) / factorial(n + m)), m * phi);
#endif
  }

  inline cmplx Rnm_Nishimura2002_reversed(const int n, const int m) const {
    const double cos_theta = -std::get<2>(X) / rho;
#if defined(_LIBCPP_VERSION)
    return std::polar(bem_sf::assoc_legendre(n, m, cos_theta) * (std::pow(rho, n) / factorial(n + m)), m * phi);
#else
    return std::polar(std::assoc_legendre(n, m, cos_theta) * (std::pow(rho, n) / factorial(n + m)), m * phi);
#endif
  }

  inline cmplx Snm_Nishimura2002(const int n, const int m) const {
    const double cos_theta = std::get<2>(X) / rho;
#if defined(_LIBCPP_VERSION)
    return std::polar(bem_sf::assoc_legendre(n, m, cos_theta) * (factorial(n - m) / std::pow(rho, n + 1)), m * phi);
#else
    return std::polar(std::assoc_legendre(n, m, cos_theta) * (factorial(n - m) / std::pow(rho, n + 1)), m * phi);
#endif
  }

  inline cmplx Snm_Nishimura2002_reversed(const int n, const int m) const {
    const double cos_theta = -std::get<2>(X) / rho;
#if defined(_LIBCPP_VERSION)
    return std::polar(bem_sf::assoc_legendre(n, m, cos_theta) * (factorial(n - m) / std::pow(rho, n + 1)), m * phi);
#else
    return std::polar(std::assoc_legendre(n, m, cos_theta) * (factorial(n - m) / std::pow(rho, n + 1)), m * phi);
#endif
  }

#define GreenGardAndRokhlin1997

public:
  /* ---------------------- P2M ------------------- */

  inline cmplx Ypq(const int n, int m) const {
    if (std::abs(m) > n || n < 0)
      return {0.0, 0.0};
    else if ((pre_mask & PrecomputeYpq) && n <= pre_max_ypq && std::abs(m) <= pre_max_ypq)
      return this->pre_compute_div_rhon1_i_absk_A_FMM[n][m + pre_max];
    else
      return this->sph_harmonics(n, m) / ((std::pow(rho, n + 1)) * static_cast<std::complex<double>>(i_absk_A_FMM[n][m + N_i_absk_A_FMM]));
  }

  inline std::array<cmplx, 2> p2mFunction(const int n, int m, const std::array<double, 3>& normal) const {
#ifdef GreenGardAndRokhlin1997
    m = -m; //! マイナスをつけるわけは，モーメントとなる側Y(a,b)_n^-mだからだ．eq. (5.3)
    cmplx Rnm = this->sph_harmonics_rho(n, m);
    auto [J0, J1, J2] = Jacobian_inv();
    const double cos_theta = std::get<2>(X) / rho;
    const int abs_m = std::abs(m);
    cmplx grad_Rnm1_dot_normal = (static_cast<double>(n) * (Rnm / rho)) * Dot(normal, J0); //! rについての微分
#if defined(_LIBCPP_VERSION)
    grad_Rnm1_dot_normal += (std::get<2>(X) / r2D * abs_m - (bem_sf::assoc_legendre(n, abs_m + 1, cos_theta) / bem_sf::assoc_legendre(n, abs_m, cos_theta))) * Rnm * Dot(normal, J1); //! thetaについての微分
#else
    grad_Rnm1_dot_normal += (std::get<2>(X) / r2D * abs_m - (std::assoc_legendre(n, abs_m + 1, cos_theta) / std::assoc_legendre(n, abs_m, cos_theta))) * Rnm * Dot(normal, J1); //! thetaについての微分
#endif
    grad_Rnm1_dot_normal += (cmplx(0., m) * Rnm) * Dot(normal, J2); //! phiについての微分
    return {Rnm, grad_Rnm1_dot_normal};
#elif defined(Nishimura2002)
    const double cos_theta = std::get<2>(X) / rho;
    cmplx Rnm = this->Rnm_Nishimura2002(n, m);
    auto [J0, J1, J2] = Jacobian_inv();
    const int abs_m = std::abs(m);
    cmplx grad_Rnm1_dot_normal = (static_cast<double>(n) * (Rnm / rho)) * Dot(normal, J0); //! rについての微分
#if defined(_LIBCPP_VERSION)
    grad_Rnm1_dot_normal += (std::get<2>(X) / r2D * abs_m - (bem_sf::assoc_legendre(n, abs_m + 1, cos_theta) / bem_sf::assoc_legendre(n, abs_m, cos_theta))) * Rnm * Dot(normal, J1); //! thetaについての微分
#else
    grad_Rnm1_dot_normal += (std::get<2>(X) / r2D * abs_m - (std::assoc_legendre(n, abs_m + 1, cos_theta) / std::assoc_legendre(n, abs_m, cos_theta))) * Rnm * Dot(normal, J1); //! thetaについての微分
#endif
    grad_Rnm1_dot_normal += (cmplx(0., m) * Rnm) * Dot(normal, J2); //! phiについての微分
    return {Rnm, grad_Rnm1_dot_normal};
#else
    static_assert(false, "Please define either GreenGardAndRokhlin1997 or Nishimura2002");
#endif
  }

  /* ---------------------- M2M ------------------- */

  inline cmplx m2mFunction(const int j, const int k, const int n, const int m) const {
#ifdef GreenGardAndRokhlin1997
    return AAA_M2M_FMM[j][k + N_AAA_M2M_FMM][n][m + N_AAA_M2M_FMM] * this->sph_harmonics_rho(n, -m);
#elif defined(Nishimura2002)
    return Rnm_Nishimura2002(n, m);
#else
    static_assert(false, "Please define either GreenGardAndRokhlin1997 or Nishimura2002");
#endif
  }

  /* ---------------------- M2L ------------------- */

  inline cmplx m2lFunction(const int j, const int k, const int n, const int m) const {
#ifdef GreenGardAndRokhlin1997
    return AAA_M2L_FMM[j][k + N_AAA_M2L_FMM][n][m + N_AAA_M2L_FMM] * this->sph_harmonics_div_rhon1(j + n, m - k);
#elif defined(Nishimura2002)
    return static_cast<double>(negOnePower(j)) * std::conj(Snm_Nishimura2002_reversed(j + n, k + m));
#else
    static_assert(false, "Please define either GreenGardAndRokhlin1997 or Nishimura2002");
#endif
  }

  // inline cmplx m2lFunction_Convolution_farPart(const int p, const int q) const {
  //    return AAA_M2L_FMM[j][k + N_AAA_M2L_FMM][n][m + N_AAA_M2L_FMM] * this->sph_harmonics_div_rhon1(j + n, m - k);
  // }

  // inline cmplx m2lFunction_Convolution_map2L(const int j, const int k) const {
  //    return AAA_M2L_FMM[j][k + N_AAA_M2L_FMM][n][m + N_AAA_M2L_FMM] * this->sph_harmonics_div_rhon1(j + n, m - k);
  // }

  /* ---------------------- L2L ------------------- */

  inline cmplx l2lFunction(const int j, const int k, const int n, const int m) const {
#ifdef GreenGardAndRokhlin1997
    return AAA_L2L_FMM[j][k + N_AAA_L2L_FMM][n][m + N_AAA_L2L_FMM] * this->sph_harmonics_rho(n - j, m - k);
#elif defined(Nishimura2002)
    return Rnm_Nishimura2002_reversed(n - j, m - k);
#else
    static_assert(false, "Please define either GreenGardAndRokhlin1997 or Nishimura2002");
#endif
  }

  /* ---------------------- L2P ------------------- */

  inline cmplx l2pFunction(const int n, const int m) const {
#ifdef GreenGardAndRokhlin1997
    return this->sph_harmonics_rho(n, m);
#elif defined(Nishimura2002)
    return this->Rnm_Nishimura2002_reversed(n, m);
#else
    static_assert(false, "Please define either GreenGardAndRokhlin1997 or Nishimura2002");
#endif
  }
};

#include <functional>

#include "basic_exception.hpp"

struct target4FMM {
  target4FMM() = default;
  target4FMM(const std::array<double, 3>& Xtarget) : Xtarget(Xtarget) {};

  std::array<double, 3> Xtarget{0., 0., 0.}; //! 位置ベクトル
  // Near-field cache (direct integration): SoA for SIMD/cache efficiency
  // Indices correspond to the global dense vector (NetworkPoint::f2Index)
  std::vector<int32_t> near_indices;
  std::vector<double> near_weights_phi;  // multiplies phin
  std::vector<double> near_weights_phin; // multiplies phi
  //
  std::vector<std::vector<int32_t>> near_indices_per_cell;     // cellごとに分けたnear_indices
  std::vector<std::vector<double>> near_weights_phi_per_cell;  // multiplies phin
  std::vector<std::vector<double>> near_weights_phin_per_cell; // multiplies phi
  // Per-cell RLE structures (for optimized access patterns within each cell)
  std::vector<std::vector<int32_t>> near_run_base_idx_per_cell;
  std::vector<std::vector<int32_t>> near_run_pos_per_cell;
  std::vector<std::vector<int32_t>> near_run_len_per_cell;
  //
  // Run-length encoding of `near_indices` to turn indirect loads into contiguous runs where possible.
  // For each run r:
  //   idx = near_run_base_idx[r] + j
  //   weight pos = near_run_pos[r] + j
  // for j=0..near_run_len[r)-1
  std::vector<int32_t> near_run_base_idx;
  std::vector<int32_t> near_run_pos;
  std::vector<int32_t> near_run_len;

  size_t near_cell_count = 0; // 訪問したセル数（統計用）

  double diagonal_coefficient = 0.;
  std::vector<std::tuple<cmplx, std::array<cmplx, 2>*>> list_Ynm_rhon1_MM_; //! 大した計算量ではないのでこだわらない

  /* ---------------------------------------- */

  /*
  this->list_Ynm_rhon1_MM_を準備
  L2Pのために，localのmomentを取得し，和を取る準備をしておく
  */

  void setL2P(const auto& Bsources) {
    auto b_deepest = Bsources.getBucketAtDeepest(Xtarget);
    if (b_deepest == nullptr)
      throw std::runtime_error("target4FMM::setL2P: b_deepest is nullptr");
    auto& Mlocal = b_deepest->MomentsLocalExpansion;
    this->list_Ynm_rhon1_MM_.assign(Mlocal.nm_set.size(), {0., nullptr});
    SphericalCoordinates S(Xtarget - Mlocal.X);
    for (const auto& [n, m] : Mlocal.nm_set) {
      this->list_Ynm_rhon1_MM_[Mlocal.index(n, m)] = {S.l2pFunction(n, m), &Mlocal.MM_[Mlocal.index(n, m)]};
    }
  }

  // old name: integrateFMM
  std::array<double, 2> integrateFarField() const {
    cmplx w_G_phin = 0., w_Gn_phi = 0.;
    for (const auto& [Ynm_rhon1, MM_] : this->list_Ynm_rhon1_MM_) {
      w_G_phin += Ynm_rhon1 * std::get<0>(*MM_);
      w_Gn_phi += Ynm_rhon1 * std::get<1>(*MM_);
    }
    return {std::real(w_G_phin), std::real(w_Gn_phi)};
  }

  /* ---------------------------------------- */

  /*

  BEMから学び始めたものにとって，まずしっかり理解しておくべきことは，source点とtarget点が異なるということ，それらが具体的に何であるかである．

  source点: 境界積分方程式を離散化した際に現れる，sumに含まれる，積分お重みも含めた１項分のことである．ここではより具体的にガウス点である．
  target点: 積分方程式を評価する点であり，境界要素法では，通常，節点の位置である．

  この違いをしっかり理解しておくことが重要である．これによって以下のような疑問も把握できる：
  例えば，途中でツリーが生成をdで打ち切った場合，target点はlevel dのバケツに入り，近傍ソース点を探すことになる．
  levelが低ければ低いほど，近傍ソース点の数は多くなる傾向にある．
  ツリー構造は，ソース点の分布に応じて適切に成長させているので，続くtarget点にとっての積分計算にとっては最適でない可能性がある．

   ---

  ８分木でツリー構造は成長していく．そのバケツの中で保存するソース点が少ない場合は，
  それ以上，ツリーをさせると非効率になってしまうので，８分割しない．

  この場合，処理に注意が必要である．

  例えば，$d$層までツリーが成長したとする．
  $d$層のバケツは，上層$<d$層それぞれの層において，M2Lでモーメントを受け取っているわけだが，
  L2Lで$d$層までモーメントを運ぶことができない．
  そのため，そのような途中で止まったバケツは，そのバケツに対して直接積分を行う必要がある．

  */

  void setDirectIntegration(const auto& Bucket_sources) {
    auto b_deepest = Bucket_sources.getBucketAtDeepest(this->Xtarget);

    thread_local DirectAccumulator acc;
    acc.reset();
    near_cell_count = 0;

    auto setNear = [&](const auto* B) {
      ++near_cell_count;
      for (const auto& source : B->data1D_vector) {
        if (!source->fill_direct_entries)
          continue;
        if (source->fill_direct_entries_nonadj && !source->isAdjacentTo(this)) {
          source->fill_direct_entries_nonadj(this, acc);
        } else {
          source->fill_direct_entries(this, acc);
        }
      }
    };

    setNear(b_deepest);

    for (const auto& b : b_deepest->buckets_near)
      setNear(b);

    auto set_parent = [&](auto b, auto& set_parent__) -> void {
      if (b != nullptr) {
        for (auto B : b->buckets_near)
          if (!B->hasChildren())
            setNear(B);
        set_parent__(b->parent, set_parent__);
      }
    };
    set_parent(b_deepest->parent, set_parent);

    // Sort touched_list (much smaller than old all_entries)
    std::sort(acc.touched_list.begin(), acc.touched_list.end());

    // Extract results into sparse arrays
    near_indices.clear();
    near_weights_phi.clear();
    near_weights_phin.clear();

    const std::size_t n_touched = acc.touched_list.size();
    near_indices.reserve(n_touched);
    near_weights_phi.reserve(n_touched);
    near_weights_phin.reserve(n_touched);

    for (int32_t idx : acc.touched_list) {
      near_indices.push_back(idx);
      near_weights_phi.push_back(acc.phi_acc[idx]);
      near_weights_phin.push_back(acc.phin_acc[idx]);
      // Reset for next target
      acc.phi_acc[idx] = 0.;
      acc.phin_acc[idx] = 0.;
      acc.touched_flags[idx] = 0;
    }
    acc.update_high_water_mark();

    const std::size_t n = near_indices.size();

    // Build contiguous runs for unified arrays (near_indices is already sorted by idx).
    near_run_base_idx.clear();
    near_run_pos.clear();
    near_run_len.clear();
    if (n > 0) {
      near_run_base_idx.reserve(n / 4 + 1);
      near_run_pos.reserve(n / 4 + 1);
      near_run_len.reserve(n / 4 + 1);
      std::size_t pos0 = 0;
      int32_t base = near_indices[0];
      int32_t prev = base;
      for (std::size_t i = 1; i < n; ++i) {
        const int32_t cur = near_indices[i];
        if (cur == prev + 1) {
          prev = cur;
          continue;
        }
        const std::size_t len = i - pos0;
        near_run_base_idx.push_back(base);
        near_run_pos.push_back(static_cast<int32_t>(pos0));
        near_run_len.push_back(static_cast<int32_t>(len));
        pos0 = i;
        base = cur;
        prev = cur;
      }
      const std::size_t len = n - pos0;
      near_run_base_idx.push_back(base);
      near_run_pos.push_back(static_cast<int32_t>(pos0));
      near_run_len.push_back(static_cast<int32_t>(len));
    }
  }

  std::array<double, 2> integrateNearField(const double* phi_by_index, const double* phin_by_index) const {
    double ret0 = 0.0;
    double ret1 = 0.0;
    if (!bemNearUseRunsEnabled()) {
      const std::size_t n = near_indices.size();
      for (std::size_t i = 0; i < n; ++i) {
        const int32_t idx = near_indices[i];
        ret0 += phin_by_index[idx] * near_weights_phi[i];
        ret1 += phi_by_index[idx] * near_weights_phin[i];
      }
      return {ret0, ret1};
    }
    const std::size_t nr = near_run_len.size();
    for (std::size_t r = 0; r < nr; ++r) {
      const int32_t base = near_run_base_idx[r];
      const int32_t pos = near_run_pos[r];
      const int32_t len = near_run_len[r];
      double local0 = 0.0;
      double local1 = 0.0;
#pragma omp simd reduction(+ : local0, local1)
      for (int32_t j = 0; j < len; ++j) {
        const int32_t idx = base + j;
        local0 += phin_by_index[idx] * near_weights_phi[static_cast<std::size_t>(pos + j)];
        local1 += phi_by_index[idx] * near_weights_phin[static_cast<std::size_t>(pos + j)];
      }
      ret0 += local0;
      ret1 += local1;
    }
    return {ret0, ret1};
  }
};

/* -------------------------------------------------------------------------- */

inline constexpr std::size_t computeSizeM2M(std::size_t N) {
  std::size_t count = 0;
  for (std::size_t j = 0; j <= N; ++j)
    for (int k = -static_cast<int>(j); k <= static_cast<int>(j); ++k)
      for (std::size_t n = 0; n <= j; ++n)
        for (int m = -static_cast<int>(n); m <= static_cast<int>(n); ++m)
          if (AAA_M2M_FMM[j][k + N_AAA_M2M_FMM][n][m + N_AAA_M2M_FMM].real() != 0.0 || AAA_M2M_FMM[j][k + N_AAA_M2M_FMM][n][m + N_AAA_M2M_FMM].imag() != 0.0)
            ++count;
  return count;
}

/* -------------------------------------------------------------------------- */
template <int N>
constexpr std::array<std::array<int, 2>, (N + 1) * (N + 1)> make_nm_set() {
  std::array<std::array<int, 2>, (N + 1) * (N + 1)> arr{};
  int ind = 0;
  for (int n = 0; n <= N; ++n)
    for (int m = -n; m <= n; ++m)
      arr[ind++] = {n, m};
  return arr;
}

template <int N>
constexpr auto make_zero_MM() {
  std::array<std::array<cmplx, 2>, (N + 1) * (N + 1)> result{};
  for (auto& pair : result)
    pair = {0.0, 0.0};
  return result;
}

inline double fma(const double& a, const double& b, const double& c, const double& d, const double& e) { return std::fma(a, b, std::fma(c, d, e)); }

#endif
