/***************************************************
 * gauss-hermite.c
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Gauss-Hermite Quadrature integration for
 * static disorder ensemble averages.
 *
 * See p. 96 in Notebook 1.
 *
 ***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>

#include "qdas.h"
#include "aux.h"
#include "gauss-hermite.h"

/* typedef */
struct gh_sdisorder_term_struct 
{
  size_t i,j; /* element (i,j) */
  double sigma; /* standard dev. */
};

typedef struct gh_sdisorder_term_struct gh_sdisorder_term;

/* Global variables */

/* sets abscissas and weights for Hermite quadrature.
   borrow from quadrule.f90 */

void gauss_hermite_set(double *xtab, double *weight, size_t order)
{
  int i;
  double sum;
  
  /* order must be between 1 and 20 */
  
  if ( order == 1 ) {
    
    xtab[0] = 0.0E+00;
    weight[0] = 1.0;
    
  } else if ( order == 2 ) {
    
    xtab[0] = - 0.707106781186547524400844362105E+00;
    xtab[1] =   0.707106781186547524400844362105E+00;
    weight[0] = 0.886226925452758013649083741671E+00;
    weight[1] = 0.886226925452758013649083741671E+00;
    
  } else if ( order == 3 ) {
    
    xtab[0] = - 0.122474487139158904909864203735E+01;
    xtab[1] =   0.0E+00;
    xtab[2] =   0.122474487139158904909864203735E+01;
    weight[0] = 0.295408975150919337883027913890E+00;
    weight[1] = 0.118163590060367735153211165556E+01;
    weight[2] = 0.295408975150919337883027913890E+00;
    
  } else if ( order == 4 ) {
    
    xtab[0] = - 0.165068012388578455588334111112E+01;
    xtab[1] = - 0.524647623275290317884060253835E+00;
    xtab[2] =   0.524647623275290317884060253835E+00;
    xtab[3] =   0.165068012388578455588334111112E+01;
    weight[0] = 0.813128354472451771430345571899E-01;
    weight[1] = 0.804914090005512836506049184481E+00;
    weight[2] = 0.804914090005512836506049184481E+00;
    weight[3] = 0.813128354472451771430345571899E-01;
    
  } else if ( order == 5 ) {
    
    xtab[0] = - 0.202018287045608563292872408814E+01;
    xtab[1] = - 0.958572464613818507112770593893E+00;
    xtab[2] =   0.0E+00;
    xtab[3] =   0.958572464613818507112770593893E+00;
    xtab[4] =   0.202018287045608563292872408814E+01;
    weight[0] = 0.199532420590459132077434585942E-01;
    weight[1] = 0.393619323152241159828495620852E+00;
    weight[2] = 0.945308720482941881225689324449E+00;
    weight[3] = 0.393619323152241159828495620852E+00;
    weight[4] = 0.199532420590459132077434585942E-01;
    
  } else if ( order == 6 ) {
    
    xtab[0] = - 0.235060497367449222283392198706E+01;
    xtab[1] = - 0.133584907401369694971489528297E+01;
    xtab[2] = - 0.436077411927616508679215948251E+00;
    xtab[3] =   0.436077411927616508679215948251E+00;
    xtab[4] =   0.133584907401369694971489528297E+01;
    xtab[5] =   0.235060497367449222283392198706E+01;
    weight[0] = 0.453000990550884564085747256463E-02;
    weight[1] = 0.157067320322856643916311563508E+00;
    weight[2] = 0.724629595224392524091914705598E+00;
    weight[3] = 0.724629595224392524091914705598E+00;
    weight[4] = 0.157067320322856643916311563508E+00;
    weight[5] = 0.453000990550884564085747256463E-02;
    
  } else if ( order == 7 ) {
    
    xtab[0] = - 0.265196135683523349244708200652E+01;
    xtab[1] = - 0.167355162876747144503180139830E+01;
    xtab[2] = - 0.816287882858964663038710959027E+00;
    xtab[3] =   0.0E+00;
    xtab[4] =   0.816287882858964663038710959027E+00;
    xtab[5] =   0.167355162876747144503180139830E+01;
    xtab[6] =   0.265196135683523349244708200652E+01;
    weight[0] = 0.971781245099519154149424255939E-03;
    weight[1] = 0.545155828191270305921785688417E-01;
    weight[2] = 0.425607252610127800520317466666E+00;
    weight[3] = 0.810264617556807326764876563813E+00;
    weight[4] = 0.425607252610127800520317466666E+00;
    weight[5] = 0.545155828191270305921785688417E-01;
    weight[6] = 0.971781245099519154149424255939E-03;
    
  } else if ( order == 8 ) {
    
    xtab[0] = - 0.293063742025724401922350270524E+01;
    xtab[1] = - 0.198165675669584292585463063977E+01;
    xtab[2] = - 0.115719371244678019472076577906E+01;
    xtab[3] = - 0.381186990207322116854718885584E+00;
    xtab[4] =   0.381186990207322116854718885584E+00;
    xtab[5] =   0.115719371244678019472076577906E+01;
    xtab[6] =   0.198165675669584292585463063977E+01;
    xtab[7] =   0.293063742025724401922350270524E+01;
    weight[0] = 0.199604072211367619206090452544E-03;
    weight[1] = 0.170779830074134754562030564364E-01;
    weight[2] = 0.207802325814891879543258620286E+00;
    weight[3] = 0.661147012558241291030415974496E+00;
    weight[4] = 0.661147012558241291030415974496E+00;
    weight[5] = 0.207802325814891879543258620286E+00;
    weight[6] = 0.170779830074134754562030564364E-01;
    weight[7] = 0.199604072211367619206090452544E-03;

  } else if ( order == 9 ) {

    xtab[0] = - 0.319099320178152760723004779538E+01;
    xtab[1] = - 0.226658058453184311180209693284E+01;
    xtab[2] = - 0.146855328921666793166701573925E+01;
    xtab[3] = - 0.723551018752837573322639864579E+00;
    xtab[4] =   0.0E+00;
    xtab[5] =   0.723551018752837573322639864579E+00;
    xtab[6] =   0.146855328921666793166701573925E+01;
    xtab[7] =   0.226658058453184311180209693284E+01;
    xtab[8] =   0.319099320178152760723004779538E+01;
    weight[0] = 0.396069772632643819045862946425E-04;
    weight[1] = 0.494362427553694721722456597763E-02;
    weight[2] = 0.884745273943765732879751147476E-01;
    weight[3] = 0.432651559002555750199812112956E+00;
    weight[4] = 0.720235215606050957124334723389E+00;
    weight[5] = 0.432651559002555750199812112956E+00;
    weight[6] = 0.884745273943765732879751147476E-01;
    weight[7] = 0.494362427553694721722456597763E-02;
    weight[8] = 0.396069772632643819045862946425E-04;

  } else if ( order == 10 ) {

    xtab[0] =  - 0.343615911883773760332672549432E+01;
    xtab[1] =  - 0.253273167423278979640896079775E+01;
    xtab[2] =  - 0.175668364929988177345140122011E+01;
    xtab[3] =  - 0.103661082978951365417749191676E+01;
    xtab[4] =  - 0.342901327223704608789165025557E+00;
    xtab[5] =    0.342901327223704608789165025557E+00;
    xtab[6] =    0.103661082978951365417749191676E+01;
    xtab[7] =    0.175668364929988177345140122011E+01;
    xtab[8] =    0.253273167423278979640896079775E+01;
    xtab[9] =   0.343615911883773760332672549432E+01;
    weight[0] =  0.764043285523262062915936785960E-05;
    weight[1] =  0.134364574678123269220156558585E-02;
    weight[2] =  0.338743944554810631361647312776E-01;
    weight[3] =  0.240138611082314686416523295006E+00;
    weight[4] =  0.610862633735325798783564990433E+00;
    weight[5] =  0.610862633735325798783564990433E+00;
    weight[6] =  0.240138611082314686416523295006E+00;
    weight[7] =  0.338743944554810631361647312776E-01;
    weight[8] =  0.134364574678123269220156558585E-02;
    weight[9] = 0.764043285523262062915936785960E-05;

  } else if ( order == 11 ) {

    xtab[0] =  - 0.366847084655958251845837146485E+01;
    xtab[1] =  - 0.278329009978165177083671870152E+01;
    xtab[2] =  - 0.202594801582575533516591283121E+01;
    xtab[3] =  - 0.132655708449493285594973473558E+01;
    xtab[4] =  - 0.656809566882099765024611575383E+00;
    xtab[5] =    0.0E+00;
    xtab[6] =    0.656809566882099765024611575383E+00;
    xtab[7] =    0.132655708449493285594973473558E+01;
    xtab[8] =    0.202594801582575533516591283121E+01;
    xtab[9] =   0.278329009978165177083671870152E+01;
    xtab[10] =   0.366847084655958251845837146485E+01;
    weight[0] =  0.143956039371425822033088366032E-05;
    weight[1] =  0.346819466323345510643413772940E-03;
    weight[2] =  0.119113954449115324503874202916E-01;
    weight[3] =  0.117227875167708503381788649308E+00;
    weight[4] =  0.429359752356125028446073598601E+00;
    weight[5] =  0.654759286914591779203940657627E+00;
    weight[6] =  0.429359752356125028446073598601E+00;
    weight[7] =  0.117227875167708503381788649308E+00;
    weight[8] =  0.119113954449115324503874202916E-01;
    weight[9] = 0.346819466323345510643413772940E-03;
    weight[10] = 0.143956039371425822033088366032E-05;

  } else if ( order == 12 ) {

    xtab[0] =  - 0.388972489786978191927164274724E+01;
    xtab[1] =  - 0.302063702512088977171067937518E+01;
    xtab[2] =  - 0.227950708050105990018772856942E+01;
    xtab[3] =  - 0.159768263515260479670966277090E+01;
    xtab[4] =  - 0.947788391240163743704578131060E+00;
    xtab[5] =  - 0.314240376254359111276611634095E+00;
    xtab[6] =    0.314240376254359111276611634095E+00;
    xtab[7] =    0.947788391240163743704578131060E+00;
    xtab[8] =    0.159768263515260479670966277090E+01;
    xtab[9] =   0.227950708050105990018772856942E+01;
    xtab[10] =   0.302063702512088977171067937518E+01;
    xtab[11] =   0.388972489786978191927164274724E+01;
    weight[0] =  0.265855168435630160602311400877E-06;
    weight[1] =  0.857368704358785865456906323153E-04;
    weight[2] =  0.390539058462906185999438432620E-02;
    weight[3] =  0.516079856158839299918734423606E-01;
    weight[4] =  0.260492310264161129233396139765E+00;
    weight[5] =  0.570135236262479578347113482275E+00;
    weight[6] =  0.570135236262479578347113482275E+00;
    weight[7] =  0.260492310264161129233396139765E+00;
    weight[8] =  0.516079856158839299918734423606E-01;
    weight[9] = 0.390539058462906185999438432620E-02;
    weight[10] = 0.857368704358785865456906323153E-04;
    weight[11] = 0.265855168435630160602311400877E-06;

  } else if ( order == 13 ) {

    xtab[0] =  - 0.410133759617863964117891508007E+01;
    xtab[1] =  - 0.324660897837240998812205115236E+01;
    xtab[2] =  - 0.251973568567823788343040913628E+01;
    xtab[3] =  - 0.185310765160151214200350644316E+01;
    xtab[4] =  - 0.122005503659074842622205526637E+01;
    xtab[5] =  - 0.605763879171060113080537108602E+00;
    xtab[6] =    0.0E+00;
    xtab[7] =    0.605763879171060113080537108602E+00;
    xtab[8] =    0.122005503659074842622205526637E+01;
    xtab[9] =   0.185310765160151214200350644316E+01;
    xtab[10] =   0.251973568567823788343040913628E+01;
    xtab[11] =   0.324660897837240998812205115236E+01;
    xtab[12] =   0.410133759617863964117891508007E+01;
    weight[0] =  0.482573185007313108834997332342E-07;
    weight[1] =  0.204303604027070731248669432937E-04;
    weight[2] =  0.120745999271938594730924899224E-02;
    weight[3] =  0.208627752961699392166033805050E-01;
    weight[4] =  0.140323320687023437762792268873E+00;
    weight[5] =  0.421616296898543221746893558568E+00;
    weight[6] =  0.604393187921161642342099068579E+00;
    weight[7] =  0.421616296898543221746893558568E+00;
    weight[8] =  0.140323320687023437762792268873E+00;
    weight[9] = 0.208627752961699392166033805050E-01;
    weight[10] = 0.120745999271938594730924899224E-02;
    weight[11] = 0.204303604027070731248669432937E-04;
    weight[12] = 0.482573185007313108834997332342E-07;

  } else if ( order == 14 ) {

    xtab[0] =  - 0.430444857047363181262129810037E+01;
    xtab[1] =  - 0.346265693360227055020891736115E+01;
    xtab[2] =  - 0.274847072498540256862499852415E+01;
    xtab[3] =  - 0.209518325850771681573497272630E+01;
    xtab[4] =  - 0.147668273114114087058350654421E+01;
    xtab[5] =  - 0.878713787329399416114679311861E+00;
    xtab[6] =  - 0.291745510672562078446113075799E+00;
    xtab[7] =    0.291745510672562078446113075799E+00;
    xtab[8] =    0.878713787329399416114679311861E+00;
    xtab[9] =   0.147668273114114087058350654421E+01;
    xtab[10] =   0.209518325850771681573497272630E+01;
    xtab[11] =   0.274847072498540256862499852415E+01;
    xtab[12] =   0.346265693360227055020891736115E+01;
    xtab[13] =   0.430444857047363181262129810037E+01;
    weight[0] =  0.862859116812515794532041783429E-08;
    weight[1] =  0.471648435501891674887688950105E-05;
    weight[2] =  0.355092613551923610483661076691E-03;
    weight[3] =  0.785005472645794431048644334608E-02;
    weight[4] =  0.685055342234652055387163312367E-01;
    weight[5] =  0.273105609064246603352569187026E+00;
    weight[6] =  0.536405909712090149794921296776E+00;
    weight[7] =  0.536405909712090149794921296776E+00;
    weight[8] =  0.273105609064246603352569187026E+00;
    weight[9] = 0.685055342234652055387163312367E-01;
    weight[10] = 0.785005472645794431048644334608E-02;
    weight[11] = 0.355092613551923610483661076691E-03;
    weight[12] = 0.471648435501891674887688950105E-05;
    weight[13] = 0.862859116812515794532041783429E-08;

  } else if ( order == 15 ) {

    xtab[0] =  - 0.449999070730939155366438053053E+01;
    xtab[1] =  - 0.366995037340445253472922383312E+01;
    xtab[2] =  - 0.296716692790560324848896036355E+01;
    xtab[3] =  - 0.232573248617385774545404479449E+01;
    xtab[4] =  - 0.171999257518648893241583152515E+01;
    xtab[5] =  - 0.113611558521092066631913490556E+01;
    xtab[6] =  - 0.565069583255575748526020337198E+00;
    xtab[7] =    0.0E+00;
    xtab[8] =    0.565069583255575748526020337198E+00;
    xtab[9] =   0.113611558521092066631913490556E+01;
    xtab[10] =   0.171999257518648893241583152515E+01;
    xtab[11] =   0.232573248617385774545404479449E+01;
    xtab[12] =   0.296716692790560324848896036355E+01;
    xtab[13] =   0.366995037340445253472922383312E+01;
    xtab[14] =   0.449999070730939155366438053053E+01;
    weight[0] =  0.152247580425351702016062666965E-08;
    weight[1] =  0.105911554771106663577520791055E-05;
    weight[2] =  0.100004441232499868127296736177E-03;
    weight[3] =  0.277806884291277589607887049229E-02;
    weight[4] =  0.307800338725460822286814158758E-01;
    weight[5] =  0.158488915795935746883839384960E+00;
    weight[6] =  0.412028687498898627025891079568E+00;
    weight[7] =  0.564100308726417532852625797340E+00;
    weight[8] =  0.412028687498898627025891079568E+00;
    weight[9] = 0.158488915795935746883839384960E+00;
    weight[10] = 0.307800338725460822286814158758E-01;
    weight[11] = 0.277806884291277589607887049229E-02;
    weight[12] = 0.100004441232499868127296736177E-03;
    weight[13] = 0.105911554771106663577520791055E-05;
    weight[14] = 0.152247580425351702016062666965E-08;

  } else if ( order == 16 ) {

    xtab[0] =  - 0.468873893930581836468849864875E+01;
    xtab[1] =  - 0.386944790486012269871942409801E+01;
    xtab[2] =  - 0.317699916197995602681399455926E+01;
    xtab[3] =  - 0.254620215784748136215932870545E+01;
    xtab[4] =  - 0.195178799091625397743465541496E+01;
    xtab[5] =  - 0.138025853919888079637208966969E+01;
    xtab[6] =  - 0.822951449144655892582454496734E+00;
    xtab[7] =  - 0.273481046138152452158280401965E+00;
    xtab[8] =    0.273481046138152452158280401965E+00;
    xtab[9] =   0.822951449144655892582454496734E+00;
    xtab[10] =   0.138025853919888079637208966969E+01;
    xtab[11] =   0.195178799091625397743465541496E+01;
    xtab[12] =   0.254620215784748136215932870545E+01;
    xtab[13] =   0.317699916197995602681399455926E+01;
    xtab[14] =   0.386944790486012269871942409801E+01;
    xtab[15] =   0.468873893930581836468849864875E+01;
    weight[0] =  0.265480747401118224470926366050E-09;
    weight[1] =  0.232098084486521065338749423185E-06;
    weight[2] =  0.271186009253788151201891432244E-04;
    weight[3] =  0.932284008624180529914277305537E-03;
    weight[4] =  0.128803115355099736834642999312E-01;
    weight[5] =  0.838100413989858294154207349001E-01;
    weight[6] =  0.280647458528533675369463335380E+00;
    weight[7] =  0.507929479016613741913517341791E+00;
    weight[8] =  0.507929479016613741913517341791E+00;
    weight[9] = 0.280647458528533675369463335380E+00;
    weight[10] = 0.838100413989858294154207349001E-01;
    weight[11] = 0.128803115355099736834642999312E-01;
    weight[12] = 0.932284008624180529914277305537E-03;
    weight[13] = 0.271186009253788151201891432244E-04;
    weight[14] = 0.232098084486521065338749423185E-06;
    weight[15] = 0.265480747401118224470926366050E-09;

  } else if ( order == 17 ) {

    xtab[0] =  - 0.487134519367440308834927655662E+01;
    xtab[1] =  - 0.406194667587547430689245559698E+01;
    xtab[2] =  - 0.337893209114149408338327069289E+01;
    xtab[3] =  - 0.275776291570388873092640349574E+01;
    xtab[4] =  - 0.217350282666662081927537907149E+01;
    xtab[5] =  - 0.161292431422123133311288254454E+01;
    xtab[6] =  - 0.106764872574345055363045773799E+01;
    xtab[7] =  - 0.531633001342654731349086553718E+00;
    xtab[8] =    0.0E+00;
    xtab[9] =   0.531633001342654731349086553718E+00;
    xtab[10] =   0.106764872574345055363045773799E+01;
    xtab[11] =   0.161292431422123133311288254454E+01;
    xtab[12] =   0.217350282666662081927537907149E+01;
    xtab[13] =   0.275776291570388873092640349574E+01;
    xtab[14] =   0.337893209114149408338327069289E+01;
    xtab[15] =   0.406194667587547430689245559698E+01;
    xtab[16] =   0.487134519367440308834927655662E+01;
    weight[0] =  0.458057893079863330580889281222E-10;
    weight[1] =  0.497707898163079405227863353715E-07;
    weight[2] =  0.711228914002130958353327376218E-05;
    weight[3] =  0.298643286697753041151336643059E-03;
    weight[4] =  0.506734995762753791170069495879E-02;
    weight[5] =  0.409200341495762798094994877854E-01;
    weight[6] =  0.172648297670097079217645196219E+00;
    weight[7] =  0.401826469470411956577635085257E+00;
    weight[8] =  0.530917937624863560331883103379E+00;
    weight[9] = 0.401826469470411956577635085257E+00;
    weight[10] = 0.172648297670097079217645196219E+00;
    weight[11] = 0.409200341495762798094994877854E-01;
    weight[12] = 0.506734995762753791170069495879E-02;
    weight[13] = 0.298643286697753041151336643059E-03;
    weight[14] = 0.711228914002130958353327376218E-05;
    weight[15] = 0.497707898163079405227863353715E-07;
    weight[16] = 0.458057893079863330580889281222E-10;

  } else if ( order == 18 ) {

    xtab[0] =  - 0.504836400887446676837203757885E+01;
    xtab[1] =  - 0.424811787356812646302342016090E+01;
    xtab[2] =  - 0.357376906848626607950067599377E+01;
    xtab[3] =  - 0.296137750553160684477863254906E+01;
    xtab[4] =  - 0.238629908916668600026459301424E+01;
    xtab[5] =  - 0.183553160426162889225383944409E+01;
    xtab[6] =  - 0.130092085838961736566626555439E+01;
    xtab[7] =  - 0.776682919267411661316659462284E+00;
    xtab[8] =  - 0.258267750519096759258116098711E+00;
    xtab[9] =   0.258267750519096759258116098711E+00;
    xtab[10] =   0.776682919267411661316659462284E+00;
    xtab[11] =   0.130092085838961736566626555439E+01;
    xtab[12] =   0.183553160426162889225383944409E+01;
    xtab[13] =   0.238629908916668600026459301424E+01;
    xtab[14] =   0.296137750553160684477863254906E+01;
    xtab[15] =   0.357376906848626607950067599377E+01;
    xtab[16] =   0.424811787356812646302342016090E+01;
    xtab[17] =   0.504836400887446676837203757885E+01;
    weight[0] =  0.782819977211589102925147471012E-11;
    weight[1] =  0.104672057957920824443559608435E-07;
    weight[2] =  0.181065448109343040959702385911E-05;
    weight[3] =  0.918112686792940352914675407371E-04;
    weight[4] =  0.188852263026841789438175325426E-02;
    weight[5] =  0.186400423875446519219315221973E-01;
    weight[6] =  0.973017476413154293308537234155E-01;
    weight[7] =  0.284807285669979578595606820713E+00;
    weight[8] =  0.483495694725455552876410522141E+00;
    weight[9] = 0.483495694725455552876410522141E+00;
    weight[10] = 0.284807285669979578595606820713E+00;
    weight[11] = 0.973017476413154293308537234155E-01;
    weight[12] = 0.186400423875446519219315221973E-01;
    weight[13] = 0.188852263026841789438175325426E-02;
    weight[14] = 0.918112686792940352914675407371E-04;
    weight[15] = 0.181065448109343040959702385911E-05;
    weight[16] = 0.104672057957920824443559608435E-07;
    weight[17] = 0.782819977211589102925147471012E-11;

  } else if ( order == 19 ) {

    xtab[0] =  - 0.522027169053748216460967142500E+01;
    xtab[1] =  - 0.442853280660377943723498532226E+01;
    xtab[2] =  - 0.376218735196402009751489394104E+01;
    xtab[3] =  - 0.315784881834760228184318034120E+01;
    xtab[4] =  - 0.259113378979454256492128084112E+01;
    xtab[5] =  - 0.204923170985061937575050838669E+01;
    xtab[6] =  - 0.152417061939353303183354859367E+01;
    xtab[7] =  - 0.101036838713431135136859873726E+01;
    xtab[8] =  - 0.503520163423888209373811765050E+00;
    xtab[9] =   0.0E+00;
    xtab[10] =   0.503520163423888209373811765050E+00;
    xtab[11] =   0.101036838713431135136859873726E+01;
    xtab[12] =   0.152417061939353303183354859367E+01;
    xtab[13] =   0.204923170985061937575050838669E+01;
    xtab[14] =   0.259113378979454256492128084112E+01;
    xtab[15] =   0.315784881834760228184318034120E+01;
    xtab[16] =   0.376218735196402009751489394104E+01;
    xtab[17] =   0.442853280660377943723498532226E+01;
    xtab[18] =   0.522027169053748216460967142500E+01;
    weight[0] =  0.132629709449851575185289154385E-11;
    weight[1] =  0.216305100986355475019693077221E-08;
    weight[2] =  0.448824314722312295179447915594E-06;
    weight[3] =  0.272091977631616257711941025214E-04;
    weight[4] =  0.670877521407181106194696282100E-03;
    weight[5] =  0.798886677772299020922211491861E-02;
    weight[6] =  0.508103869090520673569908110358E-01;
    weight[7] =  0.183632701306997074156148485766E+00;
    weight[8] =  0.391608988613030244504042313621E+00;
    weight[9] = 0.502974888276186530840731361096E+00;
    weight[10] = 0.391608988613030244504042313621E+00;
    weight[11] = 0.183632701306997074156148485766E+00;
    weight[12] = 0.508103869090520673569908110358E-01;
    weight[13] = 0.798886677772299020922211491861E-02;
    weight[14] = 0.670877521407181106194696282100E-03;
    weight[15] = 0.272091977631616257711941025214E-04;
    weight[16] = 0.448824314722312295179447915594E-06;
    weight[17] = 0.216305100986355475019693077221E-08;
    weight[18] = 0.132629709449851575185289154385E-11;

  } else if ( order == 20 ) {

    xtab[0] =  - 0.538748089001123286201690041068E+01;
    xtab[1] =  - 0.460368244955074427307767524898E+01;
    xtab[2] =  - 0.394476404011562521037562880052E+01;
    xtab[3] =  - 0.334785456738321632691492452300E+01;
    xtab[4] =  - 0.278880605842813048052503375640E+01;
    xtab[5] =  - 0.225497400208927552308233334473E+01;
    xtab[6] =  - 0.173853771211658620678086566214E+01;
    xtab[7] =  - 0.123407621539532300788581834696E+01;
    xtab[8] =  - 0.737473728545394358705605144252E+00;
    xtab[9] = - 0.245340708300901249903836530634E+00;
    xtab[10] =   0.245340708300901249903836530634E+00;
    xtab[11] =   0.737473728545394358705605144252E+00;
    xtab[12] =   0.123407621539532300788581834696E+01;
    xtab[13] =   0.173853771211658620678086566214E+01;
    xtab[14] =   0.225497400208927552308233334473E+01;
    xtab[15] =   0.278880605842813048052503375640E+01;
    xtab[16] =   0.334785456738321632691492452300E+01;
    xtab[17] =   0.394476404011562521037562880052E+01;
    xtab[18] =   0.460368244955074427307767524898E+01;
    xtab[19] =   0.538748089001123286201690041068E+01;
    weight[0] =  0.222939364553415129252250061603E-12;
    weight[1] =  0.439934099227318055362885145547E-09;
    weight[2] =  0.108606937076928169399952456345E-06;
    weight[3] =  0.780255647853206369414599199965E-05;
    weight[4] =  0.228338636016353967257145917963E-03;
    weight[5] =  0.324377334223786183218324713235E-02;
    weight[6] =  0.248105208874636108821649525589E-01;
    weight[7] =  0.109017206020023320013755033535E+00;
    weight[8] =  0.286675505362834129719659706228E+00;
    weight[9] = 0.462243669600610089650328639861E+00;
    weight[10] = 0.462243669600610089650328639861E+00;
    weight[11] = 0.286675505362834129719659706228E+00;
    weight[12] = 0.109017206020023320013755033535E+00;
    weight[13] = 0.248105208874636108821649525589E-01;
    weight[14] = 0.324377334223786183218324713235E-02;
    weight[15] = 0.228338636016353967257145917963E-03;
    weight[16] = 0.780255647853206369414599199965E-05;
    weight[17] = 0.108606937076928169399952456345E-06;
    weight[18] = 0.439934099227318055362885145547E-09;
    weight[19] = 0.222939364553415129252250061603E-12;

  } else {
    printf("The Gauss-Hermite %lu-point rule is not supported yet!\n",order);
  }

  /* now we are ready to normalize the total weight to 1.0 */
  sum=0.0;
  for(i=0;i<order;i++) {
    sum=sum+weight[i];
  }
  for(i=0;i<order;i++) {
    weight[i]=weight[i]/sum;
  }
  // DONE!
}

/* convert a disorder matrix into an gh_sdisorder_term representation */
void gh_sdisorder_mat2terms(gh_sdisorder_term *terms,gsl_matrix *m) 
{
  size_t ndim,i,j;
  size_t N;
  size_t count;

  ndim=m->size1;

  count=0;
  for(i=0;i<ndim;i++) {
    if(gsl_matrix_get(m,i,i) > DBL_EPSILON) {
      terms[count].i=i;
      terms[count].j=i;
      terms[count].sigma=gsl_matrix_get(m,i,i);
      count++;
    }
  }
  /* only count upper half off-diagonal terms */
  for(i=0;i<ndim;i++) {
    for(j=i+1;j<ndim;j++) {
      if(gsl_matrix_get(m,i,j) > DBL_EPSILON) {
	terms[count].i=i;
	terms[count].j=j;
	terms[count].sigma=gsl_matrix_get(m,i,j);
	count++;
      }
    }
  }
  // DONE
}


/* allocate memory space for gh_sdisorder_series;
   the static disorder matrix should be given as the input */
gh_sdisorder_series *gh_sdisorder_series_alloc(gsl_matrix *disorder, size_t order)
{
  size_t ndim,i;
  size_t N;
  size_t size;
  gh_sdisorder_series *ret;

  ndim=disorder->size1; /* dimension of Hamiltonian */

  /* count number of non-zero matrix elements in disorder;
     this is the dimension of random variables, i.e. integration */
  N=count_static_disorder_terms(disorder);

  /* order^N is the total number of points in the N dimension GH integral series */
  size=1;
  for(i=0;i<N;i++)
    size=size*order;

  /* alloc memory space */
  ret=(gh_sdisorder_series *)malloc(sizeof(gh_sdisorder_series));

  ret->order=order;
  ret->xtab=(double *)malloc(order*sizeof(double));
  ret->wtab=(double *)malloc(order*sizeof(double));
  /* set the abs points and weights for this order of GH rule */
  gauss_hermite_set(ret->xtab,ret->wtab,order);

  ret->npoints=size;
  ret->weights=(double *)malloc(size*sizeof(double));
  ret->deltaH=(gsl_matrix **)malloc(size*sizeof(gsl_matrix *));
  for(i=0;i<size;i++) {
    ret->deltaH[i]=gsl_matrix_alloc(ndim,ndim);
  }

  return ret;
}

void gh_sdisorder_series_free(gh_sdisorder_series *ghs)
{
  size_t i;

  for(i=0;i<ghs->npoints;i++)
    gsl_matrix_free(ghs->deltaH[i]);

  free(ghs->xtab);
  free(ghs->wtab);
  free(ghs->deltaH);
  free(ghs->weights);
  free(ghs);
}

/* sort a static disorder series by the weights */
void gh_sdisorder_series_sort_descend(gh_sdisorder_series *ghs) 
{
  size_t *ilist;
  double *weights;
  gsl_matrix **deltaH;
  int i;

  ilist=(size_t *)malloc(ghs->npoints*sizeof(size_t));
  weights=(double *)malloc(ghs->npoints*sizeof(double));
  deltaH=(gsl_matrix **)malloc(ghs->npoints*sizeof(gsl_matrix *));

  // sort the weights!
  sort_double_array_decend(ilist, ghs->weights, ghs->npoints);

  // get the sorted weights and deltaH (just pointers)
  for(i=0;i<ghs->npoints;i++) {
    weights[i]=ghs->weights[ilist[i]];
    deltaH[i]=ghs->deltaH[ilist[i]];
  }
  // now replace the arrays in ghs
  for(i=0;i<ghs->npoints;i++) {
    ghs->weights[i]=weights[i];
    ghs->deltaH[i]=deltaH[i];
  }

  free(ilist);
  free(weights);
  free(deltaH);
}

/* sets up the Gauss-Hermite sum by generating an array of weights and 
   delta Hamiltonian based on the input disorder matrix;
   the input *disorder should be a matrix, with the same dimension
   as the Hamiltonian, which contains the standard dev. of Gaussian static
   disorders.
*/
void gh_sdisorder_series_set(gh_sdisorder_series *ghs, const qdas_keys *keys)
{
  size_t ndim,i,j;
  size_t N;
  size_t size,count;
  gh_sdisorder_series *ret;
  gh_sdisorder_term *sdterms;
  gsl_matrix *wmat,*xmat;
  size_t *conf_vec; // use a vector to id configurations; this vector is a N digit gh_order gin-wei number

  ndim=keys->nsize; /* dimension of Hamiltonian */

  /* count number of non-zero matrix elements in disorder;
     this is the dimension of random variables, i.e. integration */
  N=count_static_disorder_terms(keys->disorder);

  // space for the static disorder terms
  sdterms=(gh_sdisorder_term *)malloc(N*sizeof(gh_sdisorder_term));
  conf_vec=(size_t *)malloc(N*sizeof(size_t));

  /* space for weights and x abscissas points for each sdisorder terms;
     Nterms x gh_order */
  wmat=gsl_matrix_alloc(N,ghs->order);
  xmat=gsl_matrix_alloc(N,ghs->order);

  /* gh_order^N is the total number of points in the N dimension GH integral series */
  size=1;
  for(i=0;i<N;i++)
    size=size*ghs->order;

  if(ghs->npoints != size) {
    printf("Can not set Gauss-Hermite static disorder series: npoints mismatch!\n");
    exit(EXIT_FAILURE);
  }

  /* identify static disorder terms */
  gh_sdisorder_mat2terms(sdterms,keys->disorder);
  /* construct weights and x points for each sdisorder terms */
  for(i=0;i<N;i++) {
    for(j=0;j<ghs->order;j++) {
      // wi'=wi
      gsl_matrix_set(wmat,i,j,ghs->wtab[j]);
      // xi'=sqrt(2)*sigma*xi
      gsl_matrix_set(xmat,i,j,sqrt(2.0)*sdterms[i].sigma*ghs->xtab[j]);
    }
  }

  for(i=0;i<N;i++)
    conf_vec[i]=0;
  // now generate each configuration and set the weight and deltaH
  for(count=0;count<ghs->npoints;count++) {
    size_t n;
    double weight;
    // conf_vec contains the corrrect configuration, we compute weight and deltaH
    // accordingly.
    ghs->weights[count]=1.0;
    gsl_matrix_set_zero(ghs->deltaH[count]);
    for(n=0;n<N;n++) {
      ghs->weights[count]=ghs->weights[count]*gsl_matrix_get(wmat,n,conf_vec[n]);
      // assign deltaH, find out the corresponding Hamilton element first
      i=sdterms[n].i;
      j=sdterms[n].j;
      gsl_matrix_set(ghs->deltaH[count],i,j,gsl_matrix_get(xmat,n,conf_vec[n]));
      gsl_matrix_set(ghs->deltaH[count],j,i,gsl_matrix_get(xmat,n,conf_vec[n])); // necessary when i != j
    } // done setting one-exciton part of deltaH

    // now also assign the two-exciton part
    gsl_matrix_fill_2es_from_1es(ghs->deltaH[count], keys->ntes, keys->tes_list);
    
    // update conf_vec like a  N digit gh_order gin-wei number
    conf_vec[0]++;
    for(n=0;n<N;n++)
      if(conf_vec[n]==ghs->order) {
	conf_vec[n+1]++;
	conf_vec[n]=0;
      }
  }

  /* sort the disorder series by weight in descending order */
  gh_sdisorder_series_sort_descend(ghs);

  /* DONE!! clean up */
  free(sdterms);
  free(conf_vec);
  gsl_matrix_free(wmat);
  gsl_matrix_free(xmat);
}

/* Show a Gauss-Hermite static disorder integral series */
void gh_sdisorder_series_print(gh_sdisorder_series *ghs)
{
  size_t i;

  printf("Gauss-Hermite Quadrature for static disorder integration:\n");
  printf("  --- Totally %d Points ---\n",ghs->npoints);
  printf("n, (Weight), deltaH\n");
  for(i=0;i<ghs->npoints;i++) {
    printf("% 6d  (%20.16f), deltaH=\n",i+1,ghs->weights[i]);
    gsl_matrix_print(ghs->deltaH[i]);
    printf("\n");
  }
}

#ifdef MAIN
int main()
{
  gh_sdisorder_series *ghs;
  qdas_keys keys;

  keys.nsize=4;
  keys.disorder=gsl_matrix_alloc(4,4);
  keys.ntes=1;
  keys.tes_list[0].label=3;
  keys.tes_list[0].site1=1;
  keys.tes_list[0].site2=2;

  gsl_matrix_set_zero(keys.disorder);
  gsl_matrix_set(keys.disorder,1,1,30.0);
  gsl_matrix_set(keys.disorder,2,2,50.0);
//  gsl_matrix_set(keys.disorder,1,2,20.0);
//  gsl_matrix_set(keys.disorder,2,1,20.0);

  ghs=gh_sdisorder_series_alloc(keys.disorder,15);

  gh_sdisorder_series_set(ghs,&keys);
  gh_sdisorder_series_print(ghs);

  gh_sdisorder_series_free(ghs);

  return 0;
}

#endif

/*
 * $Log$
 * Revision 1.8  2007/06/14 18:07:01  platin
 *   - add keyword "RSEED" that allows using a consistent random sequence
 *     in the MC sampling. This should help with the consistency of the
 *     results across jobs.
 *
 * Revision 1.7  2007/06/01 21:36:25  platin
 *
 *   - implement sorted weights in gauss-hermite module.
 *
 * Revision 1.6  2007/06/01 19:03:47  platin
 *
 *   - minor bug.
 *
 * Revision 1.5  2007/06/01 17:58:55  platin
 *
 *   - improved the Gaussian-Hermite module; now users can select between
 *     1-20 point rules can be used.
 *   - added keyword to select between MC and GH samplings.
 *
 * Revision 1.4  2007/03/11 17:55:28  platin
 *
 *   - support two-exciton states in Gauss-Hermite deltaH.
 *
 * Revision 1.3  2007/03/11 08:20:17  platin
 *
 *   - minor change.
 *
 * Revision 1.2  2007/03/09 18:45:02  platin
 *
 *   - normalize GH weights to 1 to keep consistent numbers in tnl-dm3pes output.
 *
 * Revision 1.1  2007/03/09 06:01:52  platin
 *
 *   - import new Gauss-Hermite Quadrature module that can be used
 *     for static disorder average.
 *
 *
 */
