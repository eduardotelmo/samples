
#include <vcl\vcl.h>
#include <stdio.h>
#include <conio.h>
#include <values.h>
#include <math.h>
#include <string.h>
#include <process.h>
#include <io.h>
#include "gifsave.h"
#include "mcube.h"

#pragma hdrstop

#include "tomo3dun.h"
//---------------------------------------------------------------------------
#pragma resource "*.dfm"
Tfrmmain *frmmain;
//---------------------------------------------------------------------------
__fastcall Tfrmmain::Tfrmmain(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------

TCanvas *tc; // eduardot

#define INFINITO -1
#define ESC 27
#define NN 1
#define ARESTA 1

int aresta = 3;

// Inicio de MASTER.H //

//////////////////////////////////////////////////////////////////////////
// VDLIB.H
//////////////////////////////////////////////////////////////////////////

signed int float_sin[] =
{0,
572  , 1144, 1715, 2286, 2856, 3425, 3993, 4560, 5126,
5690 , 6252, 6813, 7371, 7927, 8481, 9032, 9580,10126,
10668,11207,11743,12275,12803,13328,13848,14364,14876,
15383,15886,16384,16876,17364,17846,18323,18794,19260,
19720,20173,20621,21062,21497,21925,22347,22762,23170,
23571,23964,24351,24730,25101,25465,25821,26169,26509,
26841,27165,27481,27788,28087,28377,28659,28932,29196,
29451,29697,29934,30162,30381,30591,30791,30982,31163,
31335,31498,31650,31794,31927,32051,32165,32269,32364,
32448,32523,32588,32642,32687,32722,32747,32762,32767};

#define INTEGER_UNIT 3.051851e-5

class ponto
{
// Incluir sistema de coordenadas (negativas).
float ox, oy, oz;
float opx, opy, opz;
float pa, pb, pc;
float dx, dy, dz;
float kx, ky, kz;
float anga, angb, angc;
float sina, cosa, sinb, cosb, sinc, cosc;
float xr, yr, zr; // Recebem o retorno das transformacoes.
int numcor;
struct paleta
	{
	unsigned char r; // De 0 a 63.
	unsigned char g;
	unsigned char b;
	} pal[256];

public :
ponto();

void modo(int);
void pal_mem();

float fcos(int);
float fsin(int);

void ini_ori(int, int, int);
void ini_esc(float, float, float);
void ini_tra(int, int, int);
void ini_rot(int, int, int);
void ini_invrot(int _anga, int _angb, int _angc); // Rotacao inversa
void ini_pro(int, int, int, float, float, float);

inline void esc(float &, float &, float &);
inline void tra(float &, float &, float &);
inline void rot(float &, float &, float &);
inline void invrot(float &, float &, float &); // Rotacao inversa
inline void pro(float &, float &, float &);

inline float dist(float x0, float y0, float xf, float yf);
inline void discret(float &x, float &y);
inline void discret(float &x, float &y, float &z);

void plot2d(int, int, unsigned char);

inline void fplot2d(int, int, unsigned int);
inline void fplot3d(int, int, int, unsigned char);

void getpal(unsigned char, unsigned char &, unsigned char &, unsigned char &);
void putpal(unsigned char, unsigned char, unsigned char, unsigned char);
void fputpal(unsigned char);
unsigned char dither(unsigned char, unsigned char, unsigned char);

~ponto(); // Liberar memoria alocada dinamicamente.
};

//////////////////////////////////////////////////////////////////////////
// DEN.H
//////////////////////////////////////////////////////////////////////////

short map_version;		/* Version of this .den file                 */
short orig_min[3];		/* Dimensions of original data file          */
short orig_max[3];
short orig_len[3];
short extr_min[3];		/* Extracted portion of original file        */
short extr_max[3];		/*   (mins and maxes will be subset of       */
short extr_len[3];		/*    orig and lengths will be <= orig)      */
short map_min[3];		/* Dimensions of this map                    */
short map_max[3];		/*   (mins will be 0 in this program and     */
short map_len[3];		/*    lens may be != extr if warps > 0)      */
short map_warps;		/* Number of warps since extraction          */
										/*   (0 = none)                              */
int map_length;		/* Total number of densities in map          */
									/*   (= product of lens)                     */

class den
{
#ifndef MIN
#define MIN(a, b) ((a) > (b) ? (b) : (a))
#endif
/* constants for .den magic number (map_version field) */
#define	MAP_CUR_VERSION		1	/* current version */
#define	MAP_CUR_VERSION_SWAB	0x0100	/* byte-swapped current version */
#define MAX_READ_SIZE	8192	/* maximum # of bytes per read(2) call */

FILE *fp;
int fd;
unsigned char *data;
int swapbytes;

public:
den();
~den();
read_shorts(int fd, short *sbuf, int shortcount, int swap);
read_words(int fd, int *wbuf, int wordcount, int swap);
void tam(short &x, short &y, short &z);
void load_den(char *nome);
inline char get(int x,int y, int z);
};

//////////////////////////////////////////////////////////////////////////
// VDLIB
//////////////////////////////////////////////////////////////////////////

ponto::ponto()
{
numcor = 256;
ini_ori(0, 0, 0); ini_tra(0, 0, 0); ini_rot(0, 0, 0);
ini_esc(1, 1, 1);
//pal_mem();
}

ponto::~ponto()
{}

void ponto::ini_ori(int _ox, int _oy, int _oz)
{ox = _ox; oy = _oy; oz = _oz;}

void ponto::ini_esc(float _kx, float _ky, float _kz)
{kx = _kx; ky = _ky; kz = _kz;}

void ponto::ini_tra(int _dx, int _dy, int _dz)
{dx = _dx; dy = _dy; dz = _dz;}

void ponto::ini_pro(int _opx, int _opy, int _opz, float _a, float _b, float _c)
{
opx = _opx; opy = _opy; opz = _opz;
pa = _a; pb = _b; pc = _c;
}

void ponto::ini_rot(int _anga, int _angb, int _angc)
{
anga=_anga; angb=_angb; angc=_angc;
sina = fsin(anga); cosa = fcos(anga);
sinb = fsin(angb); cosb = fcos(angb);
sinc = fsin(angc); cosc = fcos(angc);
}

void ponto::ini_invrot(int _anga, int _angb, int _angc)
{
anga=360-_anga; angb=360-_angb; angc=360-_angc;
sina = fsin(anga); cosa = fcos(anga);
sinb = fsin(angb); cosb = fcos(angb);
sinc = fsin(angc); cosc = fcos(angc);
}

void ponto::esc(float &x, float &y, float &z)
{
x = x - ox; y = y - oy; z = z - oz;
x = kx * x; y = ky * y; z = kz * z;
x = x + ox; y = y + oy; z = z + oz;
}

void ponto::pro(float &x, float &y, float &z)
{
float w;

x = x - opx; y = y - opy; z = z - opz;
w = pa * x + pb * y + pc * z + 1.0;
x = x/w; y = y/w; z = z/w;
x = x + opx; y = y + opy; z = z + opz;
}

void ponto::tra(float &x, float &y, float &z)
{
x = x + dx; y = y + dy; z = z + dz;
}

inline void ponto::rot(float &x, float &y, float &z)
{
x = x - ox; y = y - oy; z = z - oz;

// x
yr = y * cosa - z * sina;
zr = y * sina + z * cosa;
y = yr;
z = zr;

// z
xr = x * cosc - y * sinc;
yr = x * sinc + y * cosc;
x = xr;
y = yr;

x = x + ox; y = y + oy; z = z + oz;
}

inline void ponto::invrot(float &x, float &y, float &z)
{
x = x - ox; y = y - oy; z = z - oz;

xr = x * cosc - y * sinc;
yr = x * sinc + y * cosc;
x = xr; y = yr;

xr = x * cosb + z * sinb;
zr = -x * sinb + z * cosb;
x = xr; z = zr;

yr = y * cosa - z * sina;
zr = y * sina + z * cosa;
y = yr; z = zr;

x = x + ox; y = y + oy; z = z + oz;
}

inline float ponto::dist(float x0, float y0, float xf, float yf)
{
return sqrt(pow(xf - x0, 2) + pow(yf - y0, 2));
}

inline void ponto::discret(float &x, float &y)
{
int i, j;
float xt, yt;
float d, mind;

mind = 10;
for (i = -1; i <= 1; i++)
	for (j = -1; j <= 1; j++)
    	{
		xt = ceil(x + i);
		yt = ceil(y + j);
        d = dist(x, y, xt, yt);
        if (d < mind)
        	{
            mind = d;
            x = xt;
            y = yt;
            }
        }
}

void ponto::fplot2d(int x, int y, unsigned int cor)
{
TColor tcor;

// Criar BMP ou DIB na memoria e jogar para o canvas
// TColor eh RGB :
tcor = (TColor)((cor << 16) | (cor << 8) | cor); // Grayscale
tc->Pen->Color = tcor;
tc->Brush->Color = tcor;

tc->Rectangle(x-1, y-1, x+1, y+1);
}

float ponto::fsin(int degree)
{
 degree %= 360;

 if (degree>=0 && degree <=90)
	 return (float_sin[degree] * INTEGER_UNIT);
 if (degree>90 && degree <=180)
	 return (float_sin[(180-degree)] * INTEGER_UNIT);
 if (degree>180 && degree <=270)
	 return -(float_sin[(degree-180)] * INTEGER_UNIT);
 else
	 return -(float_sin[(360-degree)] * INTEGER_UNIT);
}

float ponto::fcos(int degree)
{
 degree %= 360;

 if (degree>=0 && degree <=270)
	 return fsin(degree+90);
 else
	 return fsin(degree-270);
}

//////////////////////////////////////////////////////////////////////
// DEN
//////////////////////////////////////////////////////////////////////

den::den()
{}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
den::~den()
{}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
den::read_shorts(int fd, short *sbuf, int shortcount, int swap)
{
		int n, c;
		int bytecount = shortcount * 2;
		char tmp0, tmp1, tmp2, tmp3;
		char *buf = (char *)sbuf;

		while (bytecount > 0) {
	n = MIN(bytecount, MAX_READ_SIZE);
	if (read(fd, buf, n) != n)
			return(0);
	bytecount -= n;
	if (swap) { // Usado em PC.
			c = n / 8;
			n -= c * 8;
			for (; c > 0; c--) {
		tmp0 = buf[0]; buf[0] = buf[1]; buf[1] = tmp0;
		tmp1 = buf[2]; buf[2] = buf[3]; buf[3] = tmp1;
		tmp2 = buf[4]; buf[4] = buf[5]; buf[5] = tmp2;
		tmp3 = buf[6]; buf[6] = buf[7]; buf[7] = tmp3;
		buf += 8;
			}
			for (; n > 0; n -= 2) {
		tmp0 = buf[0]; buf[0] = buf[1]; buf[1] = tmp0;
		buf += 2;
			}
	} else {
			buf += n;
	}
		}
		return(1);
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
den::read_words(int fd, int *wbuf, int wordcount, int swap)
{
		int n, c;
		int bytecount = wordcount * 4;
		char tmp0, tmp1, tmp2, tmp3;
		char *buf = (char *)wbuf;

		while (bytecount > 0) {
	n = MIN(bytecount, MAX_READ_SIZE);
	if (read(fd, buf, n) != n)
			return(0);
	bytecount -= n;
	if (swap) {
			c = n / 8;
			n -= c * 8;
			for (; c > 0; c--) {
		tmp0 = buf[0]; buf[0] = buf[3]; buf[3] = tmp0;
		tmp1 = buf[1]; buf[1] = buf[2]; buf[2] = tmp1;
		tmp2 = buf[4]; buf[4] = buf[7]; buf[7] = tmp2;
		tmp3 = buf[5]; buf[5] = buf[6]; buf[6] = tmp3;
		buf += 8;
			}
			for (; n > 0; n -= 4) {
		tmp0 = buf[0]; buf[0] = buf[3]; buf[3] = tmp0;
		tmp1 = buf[1]; buf[1] = buf[2]; buf[2] = tmp1;
		buf += 4;
			}
	} else {
			buf += n;
	}
		}
		return(1);
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void den::tam(short &x, short &y, short &z)
{
x = extr_len[0];
y = extr_len[1];
z = extr_len[2];
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void den::load_den(char *nome)
{
fp = fopen(nome, "rb");
if (fp==NULL)
	return;
setvbuf(fp, NULL, _IOFBF, 32766);
fd = fileno(fp);

/* read the magic number */
if (!read_shorts(fd, &map_version, 1, 0))
	{
	fprintf(stderr, "read failed on file %s (empty file?)\n", nome);
	return;
	}
if (map_version == MAP_CUR_VERSION)
	{
	swapbytes = 0;
	}
else
if (map_version == MAP_CUR_VERSION_SWAB)
	{
	swapbytes = 1;
	}
else
	{
	fprintf(stderr, "file %s is not a density file\n", nome);
	return;
	}

if (!read_shorts(fd, orig_min, 3, swapbytes) ||
		!read_shorts(fd, orig_max, 3, swapbytes) ||
		!read_shorts(fd, orig_len, 3, swapbytes) ||
		!read_shorts(fd, extr_min, 3, swapbytes) ||
		!read_shorts(fd, extr_max, 3, swapbytes) ||
		!read_shorts(fd, extr_len, 3, swapbytes) ||
		!read_shorts(fd, map_min, 3, swapbytes) ||
		!read_shorts(fd, map_max, 3, swapbytes) ||
		!read_shorts(fd, map_len, 3, swapbytes) ||
		!read_shorts(fd, &map_warps, 1, swapbytes) ||
		!read_words(fd, &map_length, 1, swapbytes))
	{
	fprintf(stderr, "read failed on file %s (truncated file?)\n",nome);
	return;
	}

if (map_length != map_len[0]*map_len[1]*map_len[2])
	{
	fprintf(stderr, "density file %s has an inconsistent header\n",	nome);
	return;
	}
map_version = (unsigned short)(map_version/16.0);

int z;
register int x;
register int y;

// e.setmat(extr_len[2], extr_len[1], extr_len[0]); // Aten��o
e.setmat(extr_len[0], extr_len[1], extr_len[2]);
for (z = 0; z < extr_len[2]; z++)
	for (y = 0; y < extr_len[1]; y++)
		for (x = 0; x < extr_len[0]; x++)
			e.pmat(x, y, z, (unsigned char)getc(fp));

// Exibe tamanho :
frmmain->x->Caption = extr_len[0];
frmmain->y->Caption = extr_len[1];
frmmain->z->Caption = extr_len[2];
fclose(fp);
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
inline char den::get(int x,int y,int z)
{
return e.gmat(x, y, z);
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

#define CAMX 0
#define CAMY 1
#define CAMZ 2
#define CAMR 3
#define CAMD 4
#define MAXX 512
#define MAXY 380

struct vet
{
float x;
float y;
float z;
} n, l, o;

int dx = 1, dy = 1;
int desenha;
unsigned char far *pt;
char nome[255];
int trax=0, tray=0, traz=0, rotx=0, roty=0, rotz=0, light=0;
int minden=100, maxden=255;
int sx = 0, sy = 0;
int cam = CAMZ;

unsigned char telacor[MAXX][MAXY];
float telaz[MAXX][MAXY];
float telaz2[MAXX][MAXY];
unsigned long int telax[MAXX][MAXY];
// vet telan[MAXX][MAXY];
unsigned long int norma;

inline int norm(vet &a);
inline vet pe(float k, vet a);
inline float pi(vet a, vet b);
inline vet add(vet a, vet b);
inline vet sub(vet a, vet b);
inline float ilu(vet n, vet l);
inline float ilu(vet n, vet l, vet o,
				 float ka, float kd, float ke,
				 float light, float esp);
void vari();
void traxM();
void traxm();
void trayM();
void traym();
void trazM();
void trazm();
void rotxM();
void rotxm();
void rotyM();
void rotym();
void rotzM();
void rotzm();
void lightM();
void lightm();
void mindenM();
void mindenm();
void maxdenM();
void maxdenm();
void drawcam();
void drawcamx();
void drawcamy();
void drawcamz();
void xray();
void menu();
void draw();
float filtromedia2d(int, int, int);
unsigned char filtromedia3d(int, int, int, int);
void triang();
void load();
void save(char *);
void north();
void south();
void west();
void east();
void savetela();
void loadtela();
void pal();
void densidade();
void limpatela();
void clear();
void preview();
inline float bilerp(int maxx, int maxy, float xreal, float yreal);
inline float trilerp(int maxx, int maxy, int maxz,
		     		 float xreal, float yreal, float zreal);

////////////////////////////////////////////////////////////////////////////////

short lenx, leny, lenz;
ponto p;
den d;
// pcx z;


////////////////////////////////////////////////////////////////////////////////

inline float bilerp(int maxx, int maxy, float xreal, float yreal)
{
int x0 = int(xreal), y0 = int(yreal);
int x0i, y0i;
float dx = xreal - x0, dy = yreal - y0, omdx = 1.0 - dx, omdy = 1.0 - dy;

x0i = x0 + 1;
y0i = y0 + 1;

if ((x0 < 0) || (y0 < 0))
	return 0;

if ((x0i >= maxx) || (y0i >= maxy))
	return 0;

return omdx * omdy * telaz[x0][y0]  +
       omdx *   dy * telaz[x0i][y0] +
         dx * omdy * telaz[x0][y0i] +
         dx *   dy * telaz[x0i][y0i];
}

////////////////////////////////////////////////////////////////////////////////

inline float trilerp(int maxx, int maxy, int maxz,
					 float xreal, float yreal, float zreal)
{
int x0 = int(xreal), y0 = int(yreal), z0 = int(zreal);
int x0i, y0i, z0i;
float dx = xreal - x0, dy = yreal - y0, dz = zreal - z0,
      omdx = 1.0 - dx, omdy = 1.0 - dy, omdz = 1.0 - dz;

x0i = x0 + 1;
y0i = y0 + 1;
z0i = z0 + 1;

if ((x0 < 0) || (y0 < 0) || (z0 < 0))
	return 0;

if ((x0i >= maxx) || (y0i >= maxy) || (z0i >= maxz))
	return 0;

return omdx * omdy * omdz * d.get(x0, y0, z0) +
       omdx * omdy *   dz * d.get(x0, y0, z0i) +
       omdx *   dy * omdz * d.get(x0, y0i, z0) +
       omdx *   dy *   dz * d.get(x0, y0i, z0i) +
		 dx * omdy * omdz * d.get(x0i, y0, z0) +
		 dx * omdy *   dz * d.get(x0i, y0, z0i) +
		 dx *   dy * omdz * d.get(x0i, y0i, z0) +
		 dx *   dy *   dz * d.get(x0i, y0i, z0i);
}

////////////////////////////////////////////////////////////////////////////////

// main() foi quebrada em inicio() e fim()
void inicio()
{
Variant v;

tc = frmmain->imgtela->Canvas;
v = aresta;
frmmain->delta->Caption = v;
limpatela();

// Direcao da luz :
l.x = 0;
l.y = 0;
l.z = -1;

// Direcao do observador :
o.x = 0;
o.y = 0;
o.z = -1;

d.load_den("c:\\headsmal.den"); // !!!!!!!!!!!
// d.load_den("c:\\white.den");
d.tam(lenx ,leny, lenz);

desenha = 1;
vari();
// Gera filme : grava quadrado (0,0,150,150)
/*
desenha = 0;
char nome[256];
int i;

for (i = 0; i <= 360; i+=2)
	{
    Yield();
    rotx = i;
    roty = i;
    rotz = i;
    minden = 100;
    light = 0;
    draw();
	sprintf(nome, "c:\\mpeg\\fr%dv10.gif", i);
	save(nome);
    }
*/    
}

void fim()
{
d.~den();
// p.modo(3);
p.~ponto();
// v.~visual();
}

void pal()
{}

void densidade()
{}

void vari()
{
frmmain->lblrotx->Caption = rotx;
frmmain->lblroty->Caption = roty;
frmmain->lblrotz->Caption = rotz;
frmmain->lbltrax->Caption = trax;
frmmain->lbltray->Caption = tray;
frmmain->lbltraz->Caption = traz;
frmmain->lblminden->Caption = minden;
frmmain->lblmaxden->Caption = maxden;
frmmain->lbllight->Caption = light;
}

/////////////////////////////////////////

void traxM()
{
trax++;
if (trax >= lenx)
	trax = lenx - 1;
vari();
cam = CAMX;
drawcam();
}

void traxm()
{
trax--;
if (trax < 0)
	trax = 0;
vari();
cam = CAMX;
drawcam();
}

void trayM()
{
tray++;
if (tray >= leny)
	tray = leny - 1;
vari();
cam = CAMY;
drawcam();
}

void traym()
{
tray--;
if (tray < 0)
	tray = 0;
vari();
cam = CAMY;
drawcam();
}

void trazM()
{
traz++;
if (traz >= lenz)
	traz = lenz - 1;
vari();
cam = CAMZ;
drawcam();
}

void trazm()
{
traz--;
if (traz < 0)
	traz = 0;
vari();
cam = CAMZ;
drawcam();
}

void rotxM()
{
rotx+=45;
if (rotx >= 360)
	rotx = 0;
vari();
preview();
}

void rotxm()
{
rotx-=45;
if (rotx < 0)
	rotx = 0;
vari();
preview();
}

void rotyM()
{
roty+=45;
if (roty >= 360)
	roty = 0;
vari();
preview();
}

void rotym()
{
roty-=45;
if (roty < 0)
	roty = 0;
vari();
preview();
}

void rotzM()
{
rotz+=45;
if (rotz >= 360)
	rotz = 0;
vari();
preview();
}

void rotzm()
{
rotz-=45;
if (rotz < 0)
	rotz = 0;
vari();
preview();
}

void lightm()
{
light-=5;
if (light < -255)
	light = -255;
vari();
drawcam();
}

void lightM()
{
light+=5;
if (light >= 255)
	light = 255;
vari();
drawcam();
}

void mindenm()
{
minden-=5;
if (minden < 0)
	minden = 0;
vari();
preview();
}

void mindenM()
{
minden+=5;
if (minden >= maxden)
	minden = maxden;
vari();
preview();
}

void maxdenm()
{
maxden-=5;
if (maxden <= minden)
	maxden = minden;
vari();
preview();
}

void maxdenM()
{
maxden+=5;
if (maxden >= 255)
	maxden = 255;
vari();
preview();
}

void drawcam()
{
Screen->Cursor = crHourGlass;
switch (cam)
	{
	case CAMX : drawcamx(); break;
	case CAMY : drawcamy(); break;
	case CAMZ : drawcamz(); break;
	case CAMD : drawcamz(); break;
	case CAMR : drawcamz(); break;
	}
Screen->Cursor = crDefault;
}

void drawcamx()
{
register int y;
register int z;
int x;
float xr, yr, zr;
int yri, zri;
unsigned char dens;
int cor;

limpatela();

for (y = 0; y<MAXY; y++)
	for (x = 0; x<MAXX; x++)
		telacor[x][y] = 0;

p.ini_rot(rotx, roty, rotz);
p.ini_ori(lenx>>1, leny>>1, lenz>>1);
// hidemou();
x = lenx - trax - 1;
for (z = 0; z < lenz; z++)
	for (y = 0; y < leny; y++)
		{
		dens = d.get(x, y, z);
		if ((dens >= minden) && (dens <= maxden))
			{
			xr = x; yr = y; zr = z;
			p.rot(xr, yr, zr);
			yr += sx; zr += sy;
			if (yr < 0 || zr < 0 || yr >= MAXX || zr >= MAXY)
				continue;
            yri = ceil(yr);
            zri = ceil(zr);
			cor = dens;
			telacor[yri][zri] = dens;
			if (cor != 0)
				cor += light;
			if (cor < 0)
				cor = 0;
			if (cor > 255)
				cor = 255;

			p.fplot2d(yri, zri, cor);
			}
		}
}

void drawcamy()
{
register int x;
register int z;
int y;
float xr, yr, zr;
int xri, zri;
unsigned char dens;
int cor;

limpatela();
for (y = 0; y<MAXY; y++)
	for (x = 0; x<MAXX; x++)
		telacor[x][y] = 0;
p.ini_rot(rotx, roty, rotz);
p.ini_ori(lenx>>1, leny>>1, lenz>>1);

y = leny - tray - 1;
for (z = 0; z < lenz; z++)
	for (x = 0; x < lenx; x++)
		{
		dens = d.get(x, y, z);
		if ((dens >= minden) && (dens <= maxden))
			{
			xr = x; yr = y; zr = z;
			p.rot(xr, yr, zr);
			xr += sx; zr += sy;
			if (xr < 0 || zr < 0 || xr >= MAXX || zr >= MAXY)
				continue;
            xri = ceil(xr);
            zri = ceil(zr);
			cor = dens;
			telacor[xri][zri] = dens;
			if (cor != 0)
				cor += light;
			if (cor < 0)
				cor = 0;
			if (cor > 255)
				cor = 255;
			// p.fplot2d(xr, zr, v.dither(cor));
            p.fplot2d(xri, zri, cor);
			}
		}

}

void drawcamz()
{
register int x;
register int y;
int z;
float xr, yr, zr;
int xri, yri;
unsigned char dens;
int cor;

limpatela();

for (y = 0; y<MAXY; y++)
	for (x = 0; x<MAXX; x++)
		telacor[x][y] = 0;
p.ini_rot(rotx, roty, rotz);
p.ini_ori(lenx>>1, leny>>1, lenz>>1);

z = lenz - traz - 1;
for (y = 0; y < leny; y++)
	for (x = 0; x < lenx; x++)
		{
		dens = d.get(x, y, z);
		if ((dens >= minden) && (dens <= maxden))
			{
			xr = x; yr = y; zr = z;
			p.rot(xr, yr, zr);
			xr += sx; yr += sy;
			if (xr < 0 || yr < 0 || xr >= MAXX || yr >= MAXY)
				continue;
            xri = ceil(xr);
            yri = ceil(yr);
			cor = dens;
			telacor[xri][yri] = dens;
			if (cor != 0)
				cor += light;
			if (cor < 0)
				cor = 0;
			if (cor > 255)
				cor = 255;
			// p.fplot2d(xr, yr, v.dither(cor));
			p.fplot2d(xri, yri, cor);
			}
		}

}

void xray()
{
register int z;
register int x;
int y;
float xr, yr, zr;
int xri, yri;
int dens;
int cor;
float energy;
float decay = 0.1;
float normdens;

cam = CAMR;
p.ini_rot(rotx, roty, rotz);
p.ini_ori(lenx>>1, leny>>1, lenz>>1);


for (y = 0; y<MAXY; y++)
	for (x = 0; x<MAXX; x++)
		telax[x][y] = 0;

for (y = 0; y < leny; y++)
	for (x = 0; x < lenx; x++)
    	{
        energy = 10.0;
        for (z = 0; z < lenz; z++)
			{
			dens = d.get(x, y, z);
			if ((dens >= minden) && (dens <= maxden))
				{
				xr = x; yr = y; zr = z;
				p.rot(xr, yr, zr);
				xri = ceil(xr);
	            yri = ceil(yr);
				xri += sx; yri += sy;
				if (xri < 0 || yri < 0 || xri >= MAXX || yri >= MAXY)
					continue;
                telax[xri][yri] += energy * dens;
                normdens = dens / 255.0;
                energy -= decay * normdens;
                if (energy <= 0.0)
	               	break;
				// telax[xri][yri] += pow(dens, 0.4);
                }
			}
		}
limpatela();
norma = 0;
for (y = 0; y<MAXY; y++)
	for (x = 0; x<MAXX; x++)
		if (telax[x][y] > norma)
			norma = telax[x][y];

for (y = 0; y<MAXY; y++)
	for (x = 0; x<MAXX; x++)
		{
		if (telax[x][y] == 0)
			continue;
		cor = (unsigned int)(telax[x][y] * 255 / norma);
		if (cor != 0)
			cor += light;
		if (cor < 0)
			cor = 0;
		if (cor > 255)
			cor = 255;
        p.fplot2d(x, y, cor);
		}
}

void load()
{
if (strlen(nome) > 4)
	{
	d.load_den(nome);
	d.tam(lenx ,leny, lenz);
	trax=0; tray=0; traz=0; rotx=0; roty=0; rotz=0; minden=0; maxden=255; light=0;
	sx = 100; sy = 50;
	cam = CAMZ;
	}
loadtela();
vari();
// showmou();
}

void north()
{
sy-=5;
preview();
}

void south()
{
sy+=5;
preview();
}

void east()
{
sx+=5;
preview();
}

void west()
{
sx-=5;
preview();
}

void savetela()
{
register int x;
register int y;

for (y = 0; y<MAXY; y++)
	for (x = 0; x<MAXX; x++)
		telaz[x][y] = *(pt + MAXX * y + x);
}

void loadtela()
{
register int x;
register int y;

// pal();
for (y = 0; y<MAXY; y++)
	for (x = 0; x<MAXX; x++)
		p.fplot2d(x, y, (int)telaz[x][y]);
clear();
}

void limpatela()
{
tc->Brush->Color = clBlack;
tc->Rectangle(0, 0, MAXX, MAXY);
}

void clear()
{
register int x;
register int y;

for (y = 0; y<MAXY; y++)
	for (x = 0; x<MAXX; x++)
		telaz[x][y] = 0;
}

inline vet pe(float k, vet a)
{
a.x *= k;
a.y *= k;
a.z *= k;
return a;
}

inline vet add(vet a, vet b)
{
vet ret;

ret.x = a.x + b.x;
ret.y = a.y + b.y;
ret.z = a.z + b.z;
return ret;
}

inline vet sub(vet a, vet b)
{
vet ret;

ret.x = a.x - b.x;
ret.y = a.y - b.y;
ret.z = a.z - b.z;
return ret;
}

inline float pi(vet a, vet b)
{
return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline int norm(vet &a)
{
float piaa;
float norm_factor;

piaa = pi(a, a);
if (piaa == 0) // estouro de divisao
	return -1;
norm_factor = 1.0/sqrt(piaa);
a = pe(norm_factor, a);
}

inline float ilu(vet n, vet l)
{
// Iluminacao difusa :
// normaliza normal :
if (norm(n) == -1) // estouro de divisao
	return -1;
// normaliza luz :
if (norm(l) == -1) // estouro de divisao
	return -1;
return pi(n, l);
}

inline float ilu(vet n, vet l, vet o,
				 float ka, float kd, float ke,
                 float light, float esp)
{
float ia, id, ie;
float prod;
float norm_factor;
vet r;

// Iluminacao ambiente :
ia = light;

// Iluminacao difusa :
// normaliza normal :
if (norm(n) == -1) // estouro de divisao
	return -1;
// normaliza luz :
if (norm(l) == -1) // estouro de divisao
	return -1;
id = pi(n, l);

// norm_factor = ka + kd + ke;
norm_factor = ka + kd;
if (norm_factor <= 0.0)
	norm_factor = 1.0;
// return (ka * ia + kd * id + 1.0 - ke * ie) / norm_factor;
return (ka * ia + kd * id) / norm_factor;
}

void __fastcall Tfrmmain::btndrawClick(TObject *Sender)
{
draw();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button1Click(TObject *Sender)
{
rotxM();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button2Click(TObject *Sender)
{
rotxm();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button3Click(TObject *Sender)
{
rotyM();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button4Click(TObject *Sender)
{
rotym();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button6Click(TObject *Sender)
{
rotzm();	
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button8Click(TObject *Sender)
{
traxm();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button10Click(TObject *Sender)
{
traym();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button12Click(TObject *Sender)
{
trazm();	
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button5Click(TObject *Sender)
{
rotzM();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button7Click(TObject *Sender)
{
traxM();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button9Click(TObject *Sender)
{
trayM();	
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button11Click(TObject *Sender)
{
trazM();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button13Click(TObject *Sender)
{
mindenM();	
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button14Click(TObject *Sender)
{
maxdenM();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button21Click(TObject *Sender)
{
xray();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button15Click(TObject *Sender)
{
//load();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button16Click(TObject *Sender)
{
//save();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button17Click(TObject *Sender)
{
north();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button20Click(TObject *Sender)
{
south();	
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button19Click(TObject *Sender)
{
east();	
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button18Click(TObject *Sender)
{
west();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::FormClose(TObject *Sender, TCloseAction &Action)
{
fim();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::FormShow(TObject *Sender)
{
inicio();	
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button23Click(TObject *Sender)
{
mindenm();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button24Click(TObject *Sender)
{
maxdenm();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button22Click(TObject *Sender)
{
lightM();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button25Click(TObject *Sender)
{
lightm();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::imgtelaMouseMove(TObject *Sender, TShiftState Shift,
	int X, int Y)
{
frmmain->x->Caption = X;
frmmain->y->Caption = Y;
frmmain->z->Caption = telaz[X][Y];
frmmain->cor->Caption = telacor[X][Y];
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::deltaChange(TObject *Sender)
{
Variant v;

v = frmmain->delta->Caption;
aresta = v;
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button27Click(TObject *Sender)
{
preview();
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::UpDown1Click(TObject *Sender, TUDBtnType Button)
{
Variant v;

v = frmmain->delta->Caption;
if (Button == btNext)
	{
    aresta++;
	delta->Caption = v + 1;
    }
else
    if (aresta > 1)
    	{
        aresta--;
		delta->Caption = v - 1;
        }
}
//---------------------------------------------------------------------------
void __fastcall Tfrmmain::Button28Click(TObject *Sender)
{
triang();
}
//---------------------------------------------------------------------------

struct vertice
{
int x;
int y;
int z;
};

////////////////////////////////////////////////////////////

float filtromedia2d(int x, int y, int aresta)
{
int xt, yt;
int i, j;
float soma;
float z;
int cnt;

cnt = 0;
soma = 0;
for (i = -aresta; i <= aresta; i++)
	for (j = -aresta; j <= aresta; j++)
    	{
        xt = x + i;
        yt = y + j;
        if ((xt < 0) || (xt >= MAXX) || (yt < 0) || (yt >= MAXY))
       		continue;
        /*
        if (telaz[xt][yt] == INFINITO)
        	continue;
        */
        z = telaz[xt][yt];
        soma = soma + z;
        cnt++;
		}
if (cnt != 0)
	return soma /(float)cnt;
else
	return telaz[x][y];
}

////////////////////////////////////////////////////////////

unsigned char filtromedia3d(int x, int y, int z, int aresta)
{
// Nao apague
// Utilize para filtrar volumes e gerar novos volumes
int xt, yt, zt;
int i, j, k;
int soma;
int cnt;
int dens;

cnt = 0;
soma = 0;
for (i = -aresta; i <= aresta; i++)
	for (j = -aresta; j <= aresta; j++)
	    for (k = -aresta; k <= aresta; k++)
        	{
            xt = x + i;
            yt = y + j;
            zt = z + k;
            if ((xt < 0) || (xt >= lenx) || (yt < 0) || (yt >= leny) || (zt < 0) || (zt >= lenz))
            	continue;
            dens = d.get(xt, yt, zt);
/*
            if (dens == 0)
            	continue;
*/
            soma = soma + dens;
            cnt++;
            }
if (cnt == 0)
	return -1;
return ceil((float)soma/(float)cnt);
}

////////////////////////////////////////////////////////////

void triang()
{
register int z;
register int y;
int x;
int xc, yc, zc;
int midx, midy;
unsigned char dens;
unsigned int aresta;
int i, j, k;
int cnt;
verticed cubo[8];

Screen->Cursor = crHourGlass;

aresta = 2;
p.ini_rot(rotx, roty, rotz);
midx = lenx>>1;
midy = leny>>1;

limpatela();

for (x = 0; x < lenx; x+=aresta)
	{
	frmmain->prorender->Position = ceil(100.0*(float)x/(float)lenx);
	for (y = 0; y < leny; y+=aresta)
    	{
        for (z = 0; z < lenz; z+=aresta)
        	{
            cnt = 0;
			for (i = 0; i <= 1; i++)
            	{
	   			for (j = 0; j <= 1; j++)
                	{
					for (k = 0; k <= 1; k++)
                    	{
                        xc = x + k * aresta; yc = y + i * aresta; zc = z + j * aresta;
						if ((xc < lenx) && (yc < leny) && (zc < lenz))
                        	dens = d.get(xc, yc, zc);
                        else
                        	dens = 0;
						cubo[cnt].x = xc - midx;
                        cubo[cnt].y = yc - midy;
                        cubo[cnt].z = zc - lenz;
                        cubo[cnt].dens = dens;
						cnt++;
                        }
            		}
				}
            // Nao pode pegar pontos rotacionados porque
            // espera um cubo nao rotacionado :
        	inserecubo(cubo, minden, 0); // 0 - complexo / 1 - simples (nao funciona)
	       	}
        }
	}

Screen->Cursor = crDefault; // Fim do processamento, falta plotar (1 segundo).

frmmain->prorender->Position = 0;
savevrml("c:\\teste.wrl");
return;
}

/////////////////////////////////////////

int gpixel(int x, int y)
{
return telacor[x][y];
}

/////////////////////////////////////////

void save(char *caminho)
{
int q,                  /* Counter */
red[256],         /* Red component for each color */
green[256],       /* Green component for each color */
blue[256];        /* Blue component for each color */

for (q = 0; q < 256; q++)
        {
        red[q]   = q;
        green[q] = q;
        blue[q]  = q;
        }
GIF_Create(caminho, 150, 150, 256, 8);
for (q = 0; q < 256; q++)
        GIF_SetColor(q, red[q], green[q], blue[q]);
GIF_CompressImage(0, 0, -1, -1, gpixel);
GIF_Close();
}


//---------------------------------------------------------------------------

float interpola(float x0, float y0, float xf, float yf, float iso)
{
float a, b;

if ((yf - y0) == 0)
	return INFINITO;
a = (yf - y0) / (xf - x0);
if (a == 0)
	return INFINITO;
b = y0 - a * x0;
return (iso - b) / a;
}

//---------------------------------------------------------------------------

float min(float a, float b, float c)
{
if (a <= b)
	{
	if (a <= c)
    	return a;
    }
else
	if (b <= c)
    	return b;
return c;
}

//---------------------------------------------------------------------------

float max(float a, float b, float c)
{
if (a >= b)
	{
	if (a >= c)
    	return a;
    }
else
	if (b >= c)
    	return b;
return c;
}

//---------------------------------------------------------------------------

void preview()
{
dx = aresta;
dy = aresta;
draw();
dx = 1;
dy = 1;
}

//---------------------------------------------------------------------------

void draw() // Ray casting 1999 (24 jun)
{
POINT pol[5];
int xmax, ymax; // Janela de visao
register int i, j, k;
register int y;
register int x;
float xr, yr, zr;
int xri, yri;
unsigned char dens, dens2;
int cor;
vet direcao, raio, prox;
float step;
float numsteps;
float zplano;
int len;
float imed;

cam = CAMD;

// Inicializa matrizes :
for (x = 0; x<MAXX; x++)
	{
	for (y = 0; y<MAXY; y++)
    	{
		telacor[x][y] = 0;
		telaz[x][y] = INFINITO;
        }
    }

len = max(lenx, leny, lenz);
xmax = len;
ymax = len;
zplano = -0.707 * (float)len;
direcao.x = 0;
direcao.y = 0;
direcao.z = 1;
p.ini_ori(0, 0, 0);
p.ini_rot(rotx, roty, rotz);
p.rot(direcao.x, direcao.y, direcao.z);
norm(direcao);
step = max(fabs(direcao.x), fabs(direcao.y), fabs(direcao.z));
numsteps = fabs((float)(len - zplano) / step);

p.ini_ori(lenx>>1, leny>>1, lenz>>1);
p.ini_rot(rotx, roty, rotz);

for (x = 0; x < xmax; x = x + dx)
	{
	frmmain->prorender->Position = ceil(100.0 * (float)x /(float)xmax);
	for (y = 0; y < ymax; y = y + dy)
    	{
        // Acha origem do raio
        // Substituir por codigo de rotacao de texturas 3D
        xr = (float)x;
        yr = (float)y;
        zr = zplano;
 		p.rot(xr, yr, zr);

        // Lanca raio
        for (i = 0; i <= numsteps; i++)
        	{
            raio.x = xr + (float)i * direcao.x;
            raio.y = yr + (float)i * direcao.y;
            raio.z = zr + (float)i * direcao.z;

            // Verifica se o ponto eh valido
            if (((int)raio.x < 0) || ((int)raio.y < 0) || ((int)raio.z < 0))
            	continue;
            if (((int)raio.x >= lenx) || ((int)raio.y >= leny) || ((int)raio.z >= lenz))
            	continue;
    		// Pega ponto com densidade >= isosuperficie
            // dens = trilerp(lenx, leny, lenz, raio.x, raio.y, raio.z);
			dens = d.get((int)raio.x, (int)raio.y, (int)raio.z);
			if (dens < minden)
            	continue;

            // Procura ponto com densidade < isosuperficie
            // alfa = (i-1)
            prox.x = xr + ((float)i-1) * direcao.x;
            prox.y = yr + ((float)i-1) * direcao.y;
            prox.z = zr + ((float)i-1) * direcao.z;
            // Verifica se o ponto eh valido
            if (((int)prox.x < 0) || ((int)prox.y < 0) || ((int)prox.z < 0))
            	continue;
            if (((int)prox.x >= lenx) || ((int)prox.y >= leny) || ((int)prox.z >= lenz))
            	continue;
          	// dens2 = trilerp(lenx, leny, lenz, prox.x, prox.y, prox.z);
    		dens2 = d.get((int)prox.x, (int)prox.y, (int)prox.z);
            if (dens2 <= minden)
				imed = interpola(i, dens, i-1, dens2, minden);
           	// imed = i + 0.5;
            // Ponto no sistema de coordenadas inicial
	        xr = x; // direcao.x = 0
	        yr = y; // direcao.y = 0
	        zr = zplano - imed; // zr = zplano + imed * direcao.z (direcao.z = -1.0)
            xri = ceil(xr + sx);
            yri = ceil(yr + sy);

            if ((xri < 0) || (yri < 0))
            	continue;
            if ((xri >= MAXX) || (yri >= MAXY))
            	continue;
            telacor[xri][yri] = dens;
            for (j = 0; j < dx; j++)
	            for (k = 0; k < dy; k++)
                	{
					telaz[xri+j][yri+k] = zr;
					telacor[xri+j][yri+k] = dens;
                    }
            break;
            } // i
        } // y
	} // x

limpatela();

// Desenha a imagem com iluminacao :
for (x = 0; x < lenx; x = x + dx)
	for (y = 0; y < leny; y = y + dy)
		{
        if (telaz[x][y] == INFINITO)
        	continue;
		frmmain->prorender->Position = 100.0 - ceil(100.0 * (float)x /(float)lenx);
        n.x = (telaz[x+dx][y] - telaz[x][y]) / (float)dx;
		n.y = (telaz[x][y+dy] - telaz[x][y]) / (float)dy;
		n.z = -1;
		cor = 255.0 * ilu(n, l);
		// cor = ceil(255.0 * ilu(n, l, o, 1.0, 1.0, 1.0, 0.5, 3.0));
    	cor += light;

        pol[0].x = x;
        pol[0].y = y;

        pol[1].x = x + dx + 1;
        pol[1].y = y;

        pol[2].x = x + dx + 1;
        pol[2].y = y + dy + 1;

        pol[3].x = x;
        pol[3].y = y + dy;

        pol[4].x = x;
        pol[4].y = y;

	    tc->Pen->Color = (TColor)((cor << 16) | (cor << 8) | cor); // Grayscale;
        tc->Brush->Color = (TColor)((cor << 16) | (cor << 8) | cor); // Grayscale;
        tc->Polygon(pol, 4);
		}
frmmain->prorender->Position = 0;
Screen->Cursor = crHourGlass;
Screen->Cursor = crDefault; // Fim do processamento, falta plotar (1 segundo).
}
