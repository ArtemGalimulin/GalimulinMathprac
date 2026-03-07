#include <set>
#include <cmath>
#include <gmsh.h>

int main(int argc, char** argv) {
  gmsh::initialize();
  gmsh::model::add("the hollow torus");

  double lc = 0.08;
  double R = 1.4;
  double r1 = 0.3;
  double r2 = 0.55;

  int O = gmsh::model::geo::addPoint(0, 0, 0, lc);


  // POINT A
  int A = gmsh::model::geo::addPoint(R, 0, 0, lc);

  int Au1 = gmsh::model::geo::addPoint(R, r1, 0, lc);
  int Ad1 = gmsh::model::geo::addPoint(R, -r1, 0, lc);
  int Ao1 = gmsh::model::geo::addPoint(R + r1, 0, 0, lc);
  int Ai1 = gmsh::model::geo::addPoint(R - r1, 0, 0, lc);

  int Auo1 = gmsh::model::geo::addCircleArc(Ao1, A, Au1);
  int Aui1 = gmsh::model::geo::addCircleArc(Au1, A, Ai1);
  int Aid1 = gmsh::model::geo::addCircleArc(Ai1, A, Ad1);
  int Ado1 = gmsh::model::geo::addCircleArc(Ad1, A, Ao1);

  int Au2 = gmsh::model::geo::addPoint(R, r2, 0, lc);
  int Ad2 = gmsh::model::geo::addPoint(R, -r2, 0, lc);
  int Ao2 = gmsh::model::geo::addPoint(R + r2, 0, 0, lc);
  int Ai2 = gmsh::model::geo::addPoint(R - r2, 0, 0, lc);

  int Auo2 = gmsh::model::geo::addCircleArc(Ao2, A, Au2);
  int Aui2 = gmsh::model::geo::addCircleArc(Au2, A, Ai2);
  int Aid2 = gmsh::model::geo::addCircleArc(Ai2, A, Ad2);
  int Ado2 = gmsh::model::geo::addCircleArc(Ad2, A, Ao2);

  int Acir1 = gmsh::model::geo::addCurveLoop({Auo1, Aui1, Aid1, Ado1});
  int Acir2 = gmsh::model::geo::addCurveLoop({Auo2, Aui2, Aid2, Ado2});
  int Asur = gmsh::model::geo::addPlaneSurface({Acir2, -Acir1});


  // POINT B
  int B = gmsh::model::geo::addPoint(0, 0, R, lc);

  int Bu1 = gmsh::model::geo::addPoint(0, r1, R, lc);
  int Bd1 = gmsh::model::geo::addPoint(0, -r1, R, lc);
  int Bo1 = gmsh::model::geo::addPoint(0, 0, R + r1, lc);
  int Bi1 = gmsh::model::geo::addPoint(0, 0, R - r1, lc);

  int Buo1 = gmsh::model::geo::addCircleArc(Bo1, B, Bu1);
  int Bui1 = gmsh::model::geo::addCircleArc(Bu1, B, Bi1);
  int Bid1 = gmsh::model::geo::addCircleArc(Bi1, B, Bd1);
  int Bdo1 = gmsh::model::geo::addCircleArc(Bd1, B, Bo1);

  int Bu2 = gmsh::model::geo::addPoint(0, r2, R, lc);
  int Bd2 = gmsh::model::geo::addPoint(0, -r2, R, lc);
  int Bo2 = gmsh::model::geo::addPoint(0, 0, R + r2, lc);
  int Bi2 = gmsh::model::geo::addPoint(0, 0, R - r2, lc);

  int Buo2 = gmsh::model::geo::addCircleArc(Bo2, B, Bu2);
  int Bui2 = gmsh::model::geo::addCircleArc(Bu2, B, Bi2);
  int Bid2 = gmsh::model::geo::addCircleArc(Bi2, B, Bd2);
  int Bdo2 = gmsh::model::geo::addCircleArc(Bd2, B, Bo2);

  int Bcir1 = gmsh::model::geo::addCurveLoop({Buo1, Bui1, Bid1, Bdo1});
  int Bcir2 = gmsh::model::geo::addCurveLoop({Buo2, Bui2, Bid2, Bdo2});
  int Bsur = gmsh::model::geo::addPlaneSurface({Bcir2, -Bcir1});


  // POINT C
  int C = gmsh::model::geo::addPoint(-R, 0, 0, lc);

  int Cu1 = gmsh::model::geo::addPoint(-R, r1, 0, lc);
  int Cd1 = gmsh::model::geo::addPoint(-R, -r1, 0, lc);
  int Co1 = gmsh::model::geo::addPoint(-R - r1, 0, 0, lc);
  int Ci1 = gmsh::model::geo::addPoint(-R + r1, 0, 0, lc);

  int Cuo1 = gmsh::model::geo::addCircleArc(Co1, C, Cu1);
  int Cui1 = gmsh::model::geo::addCircleArc(Cu1, C, Ci1);
  int Cid1 = gmsh::model::geo::addCircleArc(Ci1, C, Cd1);
  int Cdo1 = gmsh::model::geo::addCircleArc(Cd1, C, Co1);

  int Cu2 = gmsh::model::geo::addPoint(-R, r2, 0, lc);
  int Cd2 = gmsh::model::geo::addPoint(-R, -r2, 0, lc);
  int Co2 = gmsh::model::geo::addPoint(-R - r2, 0, 0, lc);
  int Ci2 = gmsh::model::geo::addPoint(-R + r2, 0, 0, lc);

  int Cuo2 = gmsh::model::geo::addCircleArc(Co2, C, Cu2);
  int Cui2 = gmsh::model::geo::addCircleArc(Cu2, C, Ci2);
  int Cid2 = gmsh::model::geo::addCircleArc(Ci2, C, Cd2);
  int Cdo2 = gmsh::model::geo::addCircleArc(Cd2, C, Co2);

  int Ccir1 = gmsh::model::geo::addCurveLoop({Cuo1, Cui1, Cid1, Cdo1});
  int Ccir2 = gmsh::model::geo::addCurveLoop({Cuo2, Cui2, Cid2, Cdo2});
  int Csur = gmsh::model::geo::addPlaneSurface({Ccir2, -Ccir1});


  // POINT D
  int D = gmsh::model::geo::addPoint(0, 0, -R, lc);

  int Du1 = gmsh::model::geo::addPoint(0, r1, -R, lc);
  int Dd1 = gmsh::model::geo::addPoint(0, -r1, -R, lc);
  int Do1 = gmsh::model::geo::addPoint(0, 0, -R - r1, lc);
  int Di1 = gmsh::model::geo::addPoint(0, 0, -R + r1, lc);

  int Duo1 = gmsh::model::geo::addCircleArc(Do1, D, Du1);
  int Dui1 = gmsh::model::geo::addCircleArc(Du1, D, Di1);
  int Did1 = gmsh::model::geo::addCircleArc(Di1, D, Dd1);
  int Ddo1 = gmsh::model::geo::addCircleArc(Dd1, D, Do1);

  int Du2 = gmsh::model::geo::addPoint(0, r2, -R, lc);
  int Dd2 = gmsh::model::geo::addPoint(0, -r2, -R, lc);
  int Do2 = gmsh::model::geo::addPoint(0, 0, -R - r2, lc);
  int Di2 = gmsh::model::geo::addPoint(0, 0, -R + r2, lc);

  int Duo2 = gmsh::model::geo::addCircleArc(Do2, D, Du2);
  int Dui2 = gmsh::model::geo::addCircleArc(Du2, D, Di2);
  int Did2 = gmsh::model::geo::addCircleArc(Di2, D, Dd2);
  int Ddo2 = gmsh::model::geo::addCircleArc(Dd2, D, Do2);

  int Dcir1 = gmsh::model::geo::addCurveLoop({Duo1, Dui1, Did1, Ddo1});
  int Dcir2 = gmsh::model::geo::addCurveLoop({Duo2, Dui2, Did2, Ddo2});
  int Dsur = gmsh::model::geo::addPlaneSurface({Dcir2, Dcir1});


  // ARCS
  int Ou1 = gmsh::model::geo::addPoint(0, r1, 0, lc);
  int Ou2 = gmsh::model::geo::addPoint(0, r2, 0, lc);
  int Od1 = gmsh::model::geo::addPoint(0, -r1, 0, lc);
  int Od2 = gmsh::model::geo::addPoint(0, -r2, 0, lc);


  // AB
  int ABo1 = gmsh::model::geo::addCircleArc(Ao1, O, Bo1);
  int ABo2 = gmsh::model::geo::addCircleArc(Ao2, O, Bo2);
  int Ao12 = gmsh::model::geo::addLine(Ao1, Ao2);
  int Bo12 = gmsh::model::geo::addLine(Bo1, Bo2);
  int ABu1 = gmsh::model::geo::addCircleArc(Au1, Ou1, Bu1);
  int ABu2 = gmsh::model::geo::addCircleArc(Au2, Ou2, Bu2);

  int ABo = gmsh::model::geo::addCurveLoop({ABo1, Bo12, -ABo2, -Ao12});
  int ABosur = gmsh::model::geo::addPlaneSurface({ABo});

  int Au12 = gmsh::model::geo::addLine(Au1, Au2);
  int Bu12 = gmsh::model::geo::addLine(Bu1, Bu2);
  int ABu = gmsh::model::geo::addCurveLoop({ABu1, Bu12, -ABu2, -Au12});
  int ABusur = gmsh::model::geo::addSurfaceFilling({ABu});

  // BC
  int BCo1 = gmsh::model::geo::addCircleArc(Bo1, O, Co1);
  int BCo2 = gmsh::model::geo::addCircleArc(Bo2, O, Co2);
  int Co12 = gmsh::model::geo::addLine(Co1, Co2);
  int BCo = gmsh::model::geo::addCurveLoop({BCo1, Co12, -BCo2, -Bo12});
  int BCosur = gmsh::model::geo::addPlaneSurface({BCo});

  int BCu1 = gmsh::model::geo::addCircleArc(Bu1, Ou1, Cu1);
  int BCu2 = gmsh::model::geo::addCircleArc(Bu2, Ou2, Cu2);
  int Cu12 = gmsh::model::geo::addLine(Cu1, Cu2);
  int BCu = gmsh::model::geo::addCurveLoop({BCu1, Cu12, -BCu2, -Bu12});
  int BCusur = gmsh::model::geo::addSurfaceFilling({BCu});

  // CD
  int CDo1 = gmsh::model::geo::addCircleArc(Co1, O, Do1);
  int CDo2 = gmsh::model::geo::addCircleArc(Co2, O, Do2);
  int Do12 = gmsh::model::geo::addLine(Do1, Do2);
  int CDo = gmsh::model::geo::addCurveLoop({CDo1, Do12, -CDo2, -Co12});
  int CDosur = gmsh::model::geo::addPlaneSurface({CDo});

  int CDu1 = gmsh::model::geo::addCircleArc(Cu1, Ou1, Du1);
  int CDu2 = gmsh::model::geo::addCircleArc(Cu2, Ou2, Du2);
  int Du12 = gmsh::model::geo::addLine(Du1, Du2);
  int CDu = gmsh::model::geo::addCurveLoop({CDu1, Du12, -CDu2, -Cu12});
  int CDusur = gmsh::model::geo::addSurfaceFilling({CDu});

  // DA
  int DAo1 = gmsh::model::geo::addCircleArc(Do1, O, Ao1);
  int DAo2 = gmsh::model::geo::addCircleArc(Do2, O, Ao2);
  int DAo = gmsh::model::geo::addCurveLoop({DAo1, Ao12, -DAo2, -Do12});
  int DAosur = gmsh::model::geo::addPlaneSurface({DAo});

  int DAu1 = gmsh::model::geo::addCircleArc(Du1, Ou1, Au1);
  int DAu2 = gmsh::model::geo::addCircleArc(Du2, Ou2, Au2);
  int DAu = gmsh::model::geo::addCurveLoop({DAu1, Au12, -DAu2, -Du12});
  int DAusur = gmsh::model::geo::addSurfaceFilling({DAu});


  // AB down
  int ABd1 = gmsh::model::geo::addCircleArc(Ad1, Od1, Bd1);
  int ABd2 = gmsh::model::geo::addCircleArc(Ad2, Od2, Bd2);
  int Ad12 = gmsh::model::geo::addLine(Ad1, Ad2);
  int Bd12 = gmsh::model::geo::addLine(Bd1, Bd2);
  int ABd = gmsh::model::geo::addCurveLoop({ABd1, Bd12, -ABd2, -Ad12});
  int ABdsur = gmsh::model::geo::addSurfaceFilling({ABd});

  // BC down
  int BCd1 = gmsh::model::geo::addCircleArc(Bd1, Od1, Cd1);
  int BCd2 = gmsh::model::geo::addCircleArc(Bd2, Od2, Cd2);
  int Cd12 = gmsh::model::geo::addLine(Cd1, Cd2);
  int BCd = gmsh::model::geo::addCurveLoop({BCd1, Cd12, -BCd2, -Bd12});
  int BCdsur = gmsh::model::geo::addSurfaceFilling({BCd});

  // CD down
  int CDd1 = gmsh::model::geo::addCircleArc(Cd1, Od1, Dd1);
  int CDd2 = gmsh::model::geo::addCircleArc(Cd2, Od2, Dd2);
  int Dd12 = gmsh::model::geo::addLine(Dd1, Dd2);
  int CDd = gmsh::model::geo::addCurveLoop({CDd1, Dd12, -CDd2, -Cd12});
  int CDdsur = gmsh::model::geo::addSurfaceFilling({CDd});

  // DA down
  int DAd1 = gmsh::model::geo::addCircleArc(Dd1, Od1, Ad1);
  int DAd2 = gmsh::model::geo::addCircleArc(Dd2, Od2, Ad2);
  int DAd = gmsh::model::geo::addCurveLoop({DAd1, Ad12, -DAd2, -Dd12});
  int DAdsur = gmsh::model::geo::addSurfaceFilling({DAd});


  // SURFACES
  int ABi1 = gmsh::model::geo::addCircleArc(Ai1, O, Bi1);
  int ABi2 = gmsh::model::geo::addCircleArc(Ai2, O, Bi2);
  int Ai12 = gmsh::model::geo::addLine(Ai1, Ai2);
  int Bi12 = gmsh::model::geo::addLine(Bi1, Bi2);

  int BCi1 = gmsh::model::geo::addCircleArc(Bi1, O, Ci1);
  int BCi2 = gmsh::model::geo::addCircleArc(Bi2, O, Ci2);
  int Ci12 = gmsh::model::geo::addLine(Ci1, Ci2);

  int CDi1 = gmsh::model::geo::addCircleArc(Ci1, O, Di1);
  int CDi2 = gmsh::model::geo::addCircleArc(Ci2, O, Di2);
  int Di12 = gmsh::model::geo::addLine(Di1, Di2);

  int DAi1 = gmsh::model::geo::addCircleArc(Di1, O, Ai1);
  int DAi2 = gmsh::model::geo::addCircleArc(Di2, O, Ai2);


  int ABoo = gmsh::model::geo::addCurveLoop({ABo1, Bo12, -ABo2, -Ao12});
  int ABoosur = gmsh::model::geo::addSurfaceFilling({ABoo});
  int ABuu = gmsh::model::geo::addCurveLoop({ABu1, Bu12, -ABu2, -Au12});
  int ABuusur = gmsh::model::geo::addSurfaceFilling({ABuu});
  int ABdd = gmsh::model::geo::addCurveLoop({ABd1, Bd12, -ABd2, -Ad12});
  int ABddsur = gmsh::model::geo::addSurfaceFilling({ABdd});
  int ABii = gmsh::model::geo::addCurveLoop({ABi1, Bi12, -ABi2, -Ai12});
  int ABiisur = gmsh::model::geo::addSurfaceFilling({ABii});

  int ABuo = gmsh::model::geo::addCurveLoop({Auo2, ABu2, -Buo2, -ABo2});
  int ABuosur = gmsh::model::geo::addSurfaceFilling({ABuo});
  int ABui = gmsh::model::geo::addCurveLoop({Aui2, ABi2, -Bui2, -ABu2});
  int ABuisur = gmsh::model::geo::addSurfaceFilling({ABui});
  int ABid = gmsh::model::geo::addCurveLoop({Aid2, ABd2, -Bid2, -ABi2});
  int ABidsur = gmsh::model::geo::addSurfaceFilling({ABid});
  int ABdo = gmsh::model::geo::addCurveLoop({Ado2, ABo2, -Bdo2, -ABd2});
  int ABdosur = gmsh::model::geo::addSurfaceFilling({ABdo});

  int ABuoh = gmsh::model::geo::addCurveLoop({Auo1, ABu1, -Buo1, -ABo1});
  int ABuohsur = gmsh::model::geo::addSurfaceFilling({ABuoh});
  int ABuih = gmsh::model::geo::addCurveLoop({Aui1, ABi1, -Bui1, -ABu1});
  int ABuihsur = gmsh::model::geo::addSurfaceFilling({ABuih});
  int ABidh = gmsh::model::geo::addCurveLoop({Aid1, ABd1, -Bid1, -ABi1});
  int ABidhsur = gmsh::model::geo::addSurfaceFilling({ABidh});
  int ABdoh = gmsh::model::geo::addCurveLoop({Ado1, ABo1, -Bdo1, -ABd1});
  int ABdohsur = gmsh::model::geo::addSurfaceFilling({ABdoh});


  int BCoo = gmsh::model::geo::addCurveLoop({BCo1, Co12, -BCo2, -Bo12});
  int BCoosur = gmsh::model::geo::addSurfaceFilling({BCoo});
  int BCuu = gmsh::model::geo::addCurveLoop({BCu1, Cu12, -BCu2, -Bu12});
  int BCuusur = gmsh::model::geo::addSurfaceFilling({BCuu});
  int BCdd = gmsh::model::geo::addCurveLoop({BCd1, Cd12, -BCd2, -Bd12});
  int BCddsur = gmsh::model::geo::addSurfaceFilling({BCdd});
  int BCii = gmsh::model::geo::addCurveLoop({BCi1, Ci12, -BCi2, -Bi12});
  int BCiisur = gmsh::model::geo::addSurfaceFilling({BCii});

  int BCuo = gmsh::model::geo::addCurveLoop({Buo2, BCu2, -Cuo2, -BCo2});
  int BCuosur = gmsh::model::geo::addSurfaceFilling({BCuo});
  int BCui = gmsh::model::geo::addCurveLoop({Bui2, BCi2, -Cui2, -BCu2});
  int BCuisur = gmsh::model::geo::addSurfaceFilling({BCui});
  int BCid = gmsh::model::geo::addCurveLoop({Bid2, BCd2, -Cid2, -BCi2});
  int BCidsur = gmsh::model::geo::addSurfaceFilling({BCid});
  int BCdo = gmsh::model::geo::addCurveLoop({Bdo2, BCo2, -Cdo2, -BCd2});
  int BCdosur = gmsh::model::geo::addSurfaceFilling({BCdo});

  int BCuoh = gmsh::model::geo::addCurveLoop({Buo1, BCu1, -Cuo1, -BCo1});
  int BCuohsur = gmsh::model::geo::addSurfaceFilling({BCuoh});
  int BCuih = gmsh::model::geo::addCurveLoop({Bui1, BCi1, -Cui1, -BCu1});
  int BCuihsur = gmsh::model::geo::addSurfaceFilling({BCuih});
  int BCidh = gmsh::model::geo::addCurveLoop({Bid1, BCd1, -Cid1, -BCi1});
  int BCidhsur = gmsh::model::geo::addSurfaceFilling({BCidh});
  int BCdoh = gmsh::model::geo::addCurveLoop({Bdo1, BCo1, -Cdo1, -BCd1});
  int BCdohsur = gmsh::model::geo::addSurfaceFilling({BCdoh});


  int CDoo = gmsh::model::geo::addCurveLoop({CDo1, Do12, -CDo2, -Co12});
  int CDoosur = gmsh::model::geo::addSurfaceFilling({CDoo});
  int CDuu = gmsh::model::geo::addCurveLoop({CDu1, Du12, -CDu2, -Cu12});
  int CDuusur = gmsh::model::geo::addSurfaceFilling({CDuu});
  int CDdd = gmsh::model::geo::addCurveLoop({CDd1, Dd12, -CDd2, -Cd12});
  int CDddsur = gmsh::model::geo::addSurfaceFilling({CDdd});
  int CDii = gmsh::model::geo::addCurveLoop({CDi1, Di12, -CDi2, -Ci12});
  int CDiisur = gmsh::model::geo::addSurfaceFilling({CDii});

  int CDuo = gmsh::model::geo::addCurveLoop({Cuo2, CDu2, -Duo2, -CDo2});
  int CDuosur = gmsh::model::geo::addSurfaceFilling({CDuo});
  int CDui = gmsh::model::geo::addCurveLoop({Cui2, CDi2, -Dui2, -CDu2});
  int CDuisur = gmsh::model::geo::addSurfaceFilling({CDui});
  int CDid = gmsh::model::geo::addCurveLoop({Cid2, CDd2, -Did2, -CDi2});
  int CDidsur = gmsh::model::geo::addSurfaceFilling({CDid});
  int CDdo = gmsh::model::geo::addCurveLoop({Cdo2, CDo2, -Ddo2, -CDd2});
  int CDdosur = gmsh::model::geo::addSurfaceFilling({CDdo});

  int CDuoh = gmsh::model::geo::addCurveLoop({Cuo1, CDu1, -Duo1, -CDo1});
  int CDuohsur = gmsh::model::geo::addSurfaceFilling({CDuoh});
  int CDuih = gmsh::model::geo::addCurveLoop({Cui1, CDi1, -Dui1, -CDu1});
  int CDuihsur = gmsh::model::geo::addSurfaceFilling({CDuih});
  int CDidh = gmsh::model::geo::addCurveLoop({Cid1, CDd1, -Did1, -CDi1});
  int CDidhsur = gmsh::model::geo::addSurfaceFilling({CDidh});
  int CDdoh = gmsh::model::geo::addCurveLoop({Cdo1, CDo1, -Ddo1, -CDd1});
  int CDdohsur = gmsh::model::geo::addSurfaceFilling({CDdoh});


  int DAoo = gmsh::model::geo::addCurveLoop({DAo1, Ao12, -DAo2, -Do12});
  int DAoosur = gmsh::model::geo::addSurfaceFilling({DAoo});
  int DAuu = gmsh::model::geo::addCurveLoop({DAu1, Au12, -DAu2, -Du12});
  int DAuusur = gmsh::model::geo::addSurfaceFilling({DAuu});
  int DAdd = gmsh::model::geo::addCurveLoop({DAd1, Ad12, -DAd2, -Dd12});
  int DAddsur = gmsh::model::geo::addSurfaceFilling({DAdd});
  int DAii = gmsh::model::geo::addCurveLoop({DAi1, Ai12, -DAi2, -Di12});
  int DAiisur = gmsh::model::geo::addSurfaceFilling({DAii});

  int DAuo = gmsh::model::geo::addCurveLoop({Duo2, DAu2, -Auo2, -DAo2});
  int DAuosur = gmsh::model::geo::addSurfaceFilling({DAuo});
  int DAui = gmsh::model::geo::addCurveLoop({Dui2, DAi2, -Aui2, -DAu2});
  int DAuisur = gmsh::model::geo::addSurfaceFilling({DAui});
  int DAid = gmsh::model::geo::addCurveLoop({Did2, DAd2, -Aid2, -DAi2});
  int DAidsur = gmsh::model::geo::addSurfaceFilling({DAid});
  int DAdo = gmsh::model::geo::addCurveLoop({Ddo2, DAo2, -Ado2, -DAd2});
  int DAdosur = gmsh::model::geo::addSurfaceFilling({DAdo});

  int DAuoh = gmsh::model::geo::addCurveLoop({Duo1, DAu1, -Auo1, -DAo1});
  int DAuohsur = gmsh::model::geo::addSurfaceFilling({DAuoh});
  int DAuih = gmsh::model::geo::addCurveLoop({Dui1, DAi1, -Aui1, -DAu1});
  int DAuihsur = gmsh::model::geo::addSurfaceFilling({DAuih});
  int DAidh = gmsh::model::geo::addCurveLoop({Did1, DAd1, -Aid1, -DAi1});
  int DAidhsur = gmsh::model::geo::addSurfaceFilling({DAidh});
  int DAdoh = gmsh::model::geo::addCurveLoop({Ddo1, DAo1, -Ado1, -DAd1});
  int DAdohsur = gmsh::model::geo::addSurfaceFilling({DAdoh});


  // VOLUMES

  // AB
  int ABsl = gmsh::model::geo::addSurfaceLoop({
      Asur, Bsur,
      ABuosur, ABuisur, ABidsur, ABdosur,
      ABuohsur, ABuihsur, ABidhsur, ABdohsur
  });
  int ABvol = gmsh::model::geo::addVolume({ABsl});

  // BC
  int BCsl = gmsh::model::geo::addSurfaceLoop({
      Bsur, Csur,
      BCuosur, BCuisur, BCidsur, BCdosur,
      BCuohsur, BCuihsur, BCidhsur, BCdohsur
  });
  int BCvol = gmsh::model::geo::addVolume({BCsl});

  // CD
  int CDsl = gmsh::model::geo::addSurfaceLoop({
      Csur, Dsur,
      CDuosur, CDuisur, CDidsur, CDdosur,
      CDuohsur, CDuihsur, CDidhsur, CDdohsur
  });
  int CDvol = gmsh::model::geo::addVolume({CDsl});

  // DA
  int DAsl = gmsh::model::geo::addSurfaceLoop({
      Dsur, Asur,
      DAuosur, DAuisur, DAidsur, DAdosur,
      DAuohsur, DAuihsur, DAidhsur, DAdohsur
  });
  int DAvol = gmsh::model::geo::addVolume({DAsl});

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(3);
  std::set<std::string> args(argv, argv + argc);
  if (!args.count("-nopopup")) gmsh::fltk::run();
  gmsh::finalize();
  return 0;
}
