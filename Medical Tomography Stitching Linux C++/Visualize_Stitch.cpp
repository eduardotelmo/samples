//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DICOM stitching 2019
// Eduardo Telmo Fonseca Santos
// edu.telmo@gmail.com
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <bits/stdc++.h>

#define DLL_EXPORT

#include "library.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <thread>

#include <iostream>
#include <memory>
#include <thread>
#include <ctime>
#include <assert.h>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Open3D
#include <Open3D/Open3D.h>
#include <Open3D/Registration/FastGlobalRegistration.h>

// PCL
#include <pcl/io/auto_io.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/common/centroid.h>
#include <pcl/point_types.h>
#include <pcl/common/transforms.h>
#include <pcl/features/normal_3d.h>

// EIGEN
#include <Eigen/Dense>
#include <Eigen/Geometry>

// VTK-DICOM
#include <vtkDICOMMetaData.h>
#include <vtkDICOMFileSorter.h>
#include <vtkDICOMReader.h>
#include <vtkDICOMWriter.h>
#include <vtkDICOMCTGenerator.h>

// VTK
#include <vtkRenderWindow.h>
#include <vtkFileOutputWindow.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkCommand.h>
#include <vtkDICOMImageReader.h>
#include <vtkStringArray.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace open3d;
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define POINT_TYPE PointXYZRGBNormal
// #define POINT_TYPE PointXYZINormal
// #define POINT_TYPE PointXYZI

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define sizearray(a)  (sizeof(a) / sizeof((a)[0]))

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float PERC = 0.0;
float PERC_ACM = 0.0;
float PERC_OLD = 0.0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int minden = -1000; // -1000
int maxden = 32767; // 255
// int maxden = 255; // 255
double gray_compression_level = 48.0; // 1.0  (32767 grays) ... 128.0 (255 grays)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Callback for Progress Bar (graphic) - Portability
float progress_bar_graphic()
{
	float perc;

	if (PERC < 100.0)
		perc = 0.255 * registration::progress_bar_text();
	else
		perc = 0.255 * PERC;

    return perc;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Callback for Progress Bar (text) - Portability
void progress_bar_text()
{
	float perc, perc_old = -1.0;

	while(PERC < 395.0)
	{
		if (PERC < 100.0)
			perc = 0.255 * registration::progress_bar_text();
		else
			perc = 0.255 * PERC;

		if (perc != perc_old)
			{
			for(int i = 0; i < perc; i++)
				std::cout << "+";
			std::cout << endl;
			}		

		perc_old = perc;

		// usleep(100000);
		// sleep(1);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

constexpr inline double cubicInterpolate(double p[4], double x) {
	return p[1] + 0.5 * x * (p[2] - p[0] + x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] + x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
}

constexpr inline double bicubicInterpolate(double p[4][4], double x, double y) {
	double arr[4] = { 0,0,0,0 };
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);
}

constexpr inline double tricubicInterpolate(double p[4][4][4], double x, double y, double z) {
	double arr[4] = { 0,0,0,0 };
	arr[0] = bicubicInterpolate(p[0], y, z);
	arr[1] = bicubicInterpolate(p[1], y, z);
	arr[2] = bicubicInterpolate(p[2], y, z);
	arr[3] = bicubicInterpolate(p[3], y, z);
	return cubicInterpolate(arr, x);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Stitcher
{
	// Visualiza etapas em 3D
	bool visualize; //  true

	// Limiar para extrair apenas os ossos/dentes do arquivo DICOM
	unsigned int thr_bone; // 1376 

	// Decimation 
	unsigned int decimation; // 4 2

	// Escalas, dimensoes e conteudos dos arquivos DICOM [S]ource e [T]arget
	double xScaleORG, yScaleORG, zScaleORG, ScaleORG;
	double xScale, yScale, zScale, Scale;
	double xScaleS, yScaleS, zScaleS, ScaleS;
	double xScaleT, yScaleT, zScaleT, ScaleT;

	// Dimensoes das tomografias
	int* dims;
	int dimsS[3];
	int dimsT[3];
	int dimsO[3];

	// Ponteiros para os voxels das tomografias
	vtkSmartPointer<vtkImageData> sliceData;
	vtkSmartPointer<vtkImageData> sliceDataS;
	vtkSmartPointer<vtkImageData> sliceDataT;
	vtkSmartPointer<vtkImageData> sliceDataO;

	// Matriz homgenea de transformacao afim e sua matriz inversa (mudanca de coordenadas)
	Eigen::Matrix4d W;
	Eigen::Matrix4d WI;

	// Filename of INI file
	char inifile[256];

	// RANSAC parameter
	unsigned long max_iterations;
	unsigned long match_points;

	public:

	Stitcher()
	{
		// Read INI file
		strcpy(inifile, "config.ini");
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int stitch(char* file_in_1, char* file_in_2, char* file_out)
	{
		using namespace open3d;

		// utility::SetVerbosityLevel(utility::VerbosityLevel::Debug); // Exibe mensagens de DEBUG do Open3D

		std::shared_ptr<geometry::PointCloud> source, target;
		std::shared_ptr<geometry::PointCloud> sourcepp, targetpp;
		std::shared_ptr<registration::Feature> source_fpfh, target_fpfh;

		source = readDicomToOpen3D(file_in_1, 0, 255, 0);

		// Armazena variaveis relativas a [S]ource:
		dimsS[0] = dims[0]; dimsS[1] = dims[1]; dimsS[2] = dims[2];
		xScaleS = xScale; yScaleS = yScale; zScaleS = zScale; ScaleS = Scale;
		sliceDataS = sliceData;

		target = readDicomToOpen3D(file_in_2, 255, 0, 0);

		// Armazena variaveis relativas a [T]arget:
		dimsT[0] = dims[0]; dimsT[1] = dims[1]; dimsT[2] = dims[2];
		xScaleT = xScale; yScaleT = yScale; zScaleT = zScale; ScaleT = Scale;
		sliceDataT = sliceData;

		// Estima normais para visualizacao (serao calculadas de novo no preprocess do stitching)
		cout << "Estimating normals\n";
		source->EstimateNormals(open3d::geometry::KDTreeSearchParamHybrid(4, 30)); // @@@
		target->EstimateNormals(open3d::geometry::KDTreeSearchParamHybrid(4, 30)); // @@@

		// Exibe source e target
		if (visualize)
		{
			VisualizeComponent(*source);
			VisualizeComponent(*target);
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Preprocessing pointcloud (Open3D)
		cout << "\nPreprocessing\n";

		// Para testar com arquivos DICOM
		// Pre-processamento de pointcloud (Open3D): preprocess a partir da memoria
		std::tie(sourcepp, source_fpfh) = PreprocessPointCloud(source);
		std::tie(targetpp, target_fpfh) = PreprocessPointCloud(target);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Visualiza pointcloud pre-processado

		if (visualize)
		{
			cout << "\nType q to proceed\n";
			VisualizeBeforeRegistration(*sourcepp, *targetpp);
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Registration (stitching)
		// Open3D: Correspondence checker
		// Verifica correspondencia entre os pontos para eliminar parte deles

		cout << "\nCorrespondence checker\n";
		std::vector<std::reference_wrapper<const registration::CorrespondenceChecker>> correspondence_checker;

		// Todos os 3 checkers devem estar ativos e configurados em seus parametros (aqui e depois no push_back):
		auto correspondence_checker_edge_length = registration::CorrespondenceCheckerBasedOnEdgeLength(0.5); //  0.5 0.9
		auto correspondence_checker_distance = registration::CorrespondenceCheckerBasedOnDistance(11); // 11 4 2
		auto correspondence_checker_normal = registration::CorrespondenceCheckerBasedOnNormal(0.6); //  0.6 0.52359878 0.52359878

		// Pointcloud (PCL): parametros originais
		/*
		std::vector<std::reference_wrapper<const registration::CorrespondenceChecker>> correspondence_checker;
		auto correspondence_checker_edge_length = registration::CorrespondenceCheckerBasedOnEdgeLength(0.9);
		auto correspondence_checker_distance = registration::CorrespondenceCheckerBasedOnDistance(0.075);
		auto correspondence_checker_normal = registration::CorrespondenceCheckerBasedOnNormal(0.52359878);
		*/

		//  Todos os 3 checkers devem estar ativos e configurados:
		correspondence_checker.push_back(correspondence_checker_edge_length);
		correspondence_checker.push_back(correspondence_checker_distance);
		correspondence_checker.push_back(correspondence_checker_normal);

		// RANSAC - Otimizacao dos parametros da matriz de transformacao homogenea
        cout << endl << "==========================================================================" << endl << endl;
		cout << endl << ">>> RANSAC" << endl << endl;
        
	    PERC = 0.0; // !!!

		auto registration_result = registration::RegistrationRANSACBasedOnFeatureMatchingProgress(*sourcepp, *targetpp, *source_fpfh, *target_fpfh, 1.0, // 1.0 // @@@
			registration::TransformationEstimationPointToPoint(false), 4, // 4
			correspondence_checker,
			registration::RANSACConvergenceCriteria(max_iterations, match_points)); // O ultimo argumento e o numero maximo de validation points

		/////////////////////////////////
		// Quality control
		/////////////////////////////////

		cout << endl << ">>> Fitness: " << 100.0 * registration_result.fitness_ << " %" << endl;
		cout << endl << ">>> Stitching matching: " << registration::progress_total_validation() << " %" << endl;
		cout << endl << ">>> Stitching matching (raw): " << registration::total_validation() << endl;

		printf("RANSAC: Fitness %f, RMSE %f", registration_result.fitness_, registration_result.inlier_rmse_);

		float perc_match = 0.0;
		if (registration::progress_total_validation() > 0)
			{
			perc_match = 100.0 * (float)registration_result.fitness_ / (float)registration::progress_total_validation();
			cout << endl << ">>> Matching percentage: " << 100.0 * perc_match << endl;
			}
		else
			perc_match = 0.0;
	
		if (perc_match >= 0.07)
			{
			cout << endl << ">>> Stitching SUCCESS: " << perc_match << " %" << endl;
			}		
		else
			{
			cout << endl << ">>> Stitching FAILURE: " << perc_match << " %" << endl;
			}		
		
		/////////////////////////////////

		cout << endl << ">>> Progress: " << PERC << endl << endl;
		PERC_ACM = PERC;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Exibe matrizes homogeneas de transformacao afim
		W = registration_result.transformation_;
		WI = W.inverse();
		
		cout << "Transform matrix (homogeneous 4x4):\n" << W << endl;
		cout << "Inverse transform matrix (homogeneous 4x4):\n" << WI << endl;
		Eigen::Matrix4d WWI = W*WI;
		cout << "Product matrix (~Identity):\n" << WWI << endl;

		if (visualize)
		{
			cout << "\nType q to proceed\n";
			VisualizeAfterRegistration(*sourcepp, *targetpp, registration_result.transformation_);
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Transformacao dos volumes 3D para composicao (stitching):
		// 1-Calcular dimensoes do novo volume composto a partir dos volumes originais
		// 2-Calcular deslocamento dos indices para evitar indices negativos ou fora da faixa
		// 3-Fazer rotacao inversa para nao perder voxels (com interpolacao opcional)
		// 4-Fazer blending considerando que as regioes limitrofes sao menos confiaveis
		// 5-Salvar arquivo DICOM

	    cout << "TRANSFORM BEGINS" << endl;
		int dimsO[3];
		dimsO[0] = max(dimsS[0], dimsT[0]);
		dimsO[1] = max(dimsS[1], dimsT[1]);
		dimsO[2] = max(dimsS[2], dimsT[2]);

		cout << "\nSource dimensions (without transform): " << dimsS[0] << ", " << dimsS[1] << ", " << dimsS[2] << endl;
		cout << "\nOutput dimensions (without transform): " << dimsO[0] << ", " << dimsO[1] << ", " << dimsO[2] << endl;

		//  Transforma coordenadas limites do volume [S]ource para as novas posicoes (8 vertices extremos)
		double minx = DBL_MAX, miny = DBL_MAX, minz = DBL_MAX;
		double maxx = -DBL_MAX, maxy = -DBL_MAX, maxz = -DBL_MAX;
		double xS, yS, zS, xT, yT, zT;

		double minzT = DBL_MAX, maxzT = -DBL_MAX, minzSo = DBL_MAX, maxzSo = -DBL_MAX;
			
		cout << "\nMaps old Source volume coordinates into transformed coordinates (volume has 8 vertexes)\n";
		for (int i = 0; i <= 1; i++)
		{
			for (int j = 0; j <= 1; j++)
			{		
				for (int k = 0; k <= 1; k++)
				{
					if (i == 0) {
						xS = 0.0;
						xT = 0.0;
					}
					else {
						xS = dimsS[0] - 1;
						xT = dimsT[0] - 1;
					}
					if (j == 0) {
						yS = 0.0;
						yT = 0.0;
					}
					else {
						yS = dimsS[1] - 1;
						yT = dimsT[1] - 1;
					}
					if (k == 0) {
						zS = 0.0;
						zT = 0.0;
					}
					else {
						zS = dimsS[2] - 1;
						zT = dimsT[2] - 1;
					}

					Eigen::Vector4d coord_old = Eigen::Vector4d(xS, yS, zS, 1.0);
					Eigen::Vector4d coord_new = W * coord_old;
					minx = min(coord_new(0), minx);
					maxx = max(coord_new(0), maxx);
					miny = min(coord_new(1), miny);
					maxy = max(coord_new(1), maxy);
					minz = min(coord_new(2), minz);
					maxz = max(coord_new(2), maxz);

					Eigen::Vector4d coord_t = Eigen::Vector4d(xT, yT, zT, 1.0);
					minx = min(coord_t(0), minx);
					maxx = max(coord_t(0), maxx);
					miny = min(coord_t(1), miny);
					maxy = max(coord_t(1), maxy);
					minz = min(coord_t(2), minz);
					maxz = max(coord_t(2), maxz);

					minzT = min(coord_t(2), minzT);
					maxzT = max(coord_t(2), maxzT);
					minzSo = min(coord_new(2), minzSo);
					maxzSo = max(coord_new(2), maxzSo);

                    /*
					cout << "old: " << coord_old[0] << " " << coord_old[1] << " " << coord_old[2] << endl << "new: " << coord_new[0] << " " << coord_new[1] << " " << coord_new[2] << endl << endl;
                    */
				}
			}
		}

		// Deslocar os indices, subtraindo deste valor, para que o menor indice seja zero
		int x0, y0, z0;
		x0 = min(minx, 0.0);
		y0 = min(miny, 0.0);
		z0 = min(minz, 0.0);

		maxx = max(maxx, (double)dimsT[0]);
		maxy = max(maxy, (double)dimsT[1]);
		maxz = max(maxz, (double)dimsT[2]);

		dimsO[0] = ceil(maxx - minx);
		dimsO[1] = ceil(maxy - miny);
		dimsO[2] = ceil(maxz - minz);

        /*
		cout << "minx: " << minx << " miny: " << miny << " minz: " << minz << endl << " maxx: " << maxx << " maxy: " << maxy << " maxz: " << maxz << endl << endl;
		cout << "minzT: " << minzT << " maxzT: " << maxzT << " minzSo: " << minzSo << " maxzSo: " << maxzSo << endl << endl;
        */
        
		int widthT = dimsT[0];
		int heightT = dimsT[1];
		int depthT = dimsT[2];
		std::vector<uint16_t> raw_dataT(widthT* heightT* depthT, 0);

		int widthS = dimsS[0];
		int heightS = dimsS[1];
		int depthS = dimsS[2];
		std::vector<uint16_t> raw_dataS(widthS* heightS* depthS, 0);

		for (int z = 0; z < depthT; z++)
		{
			PERC = floor(100.0 * z / (float)(depthT - 1)); // !!!
			PERC = PERC + PERC_ACM;

			for (int y = 0; y < heightT; y++)
			{
				for (int x = 0; x < widthT; x++)
				{
					raw_dataT[z * heightT * widthT + y * widthT + x] = (uint16_t)std::lround(sliceDataT->GetScalarComponentAsDouble(x, y, z, 0) - minden); // !!!
				}
			}
		}

		cout << endl << ">>> Progress: " << PERC << endl << endl;
		PERC_ACM = PERC;

		for (int z = 0; z < depthS; z++)
		{
            		PERC = floor(100.0 * z / (float)(depthS - 1)); // !!!
            		PERC = PERC + PERC_ACM;

			for (int y = 0; y < heightS; y++)
			{
				for (int x = 0; x < widthS; x++)
				{
					raw_dataS[z * heightS * widthS + y * widthS + x] = (uint16_t)std::lround(sliceDataS->GetScalarComponentAsDouble(x, y, z, 0) - minden); // !!!
				}
			}
		}

		cout << endl << ">>> Progress: " << PERC << endl << endl;
        PERC_ACM = PERC;

		cout << "\nOutput dimensions (with transform): " << dimsO[0] << ", " << dimsO[1] << ", " << dimsO[2] << endl;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Menor valor de largura de borda (borda e a regiao menos confiavel proxima aos limites da tomografia): 1 voxel. 
		// Caso seja 0, desativa borda de desconfianca (considera todos os valores como igualmente confiaveis).
		unsigned int border = 10; // 10 (suavizou as emendas) // 10 // 1 // 0

		sliceDataO = vtkSmartPointer<vtkImageData>::New(); // [O]utput
		auto img = vtkSmartPointer <vtkImageData>::New();

		// Assumi os fatores de escala do [S]ource
		double spacingDataO[3];

		// Metadata
		// Create a generator for CT images.
		vtkNew<vtkDICOMCTGenerator> generator; // !!! NOVO
		// Create a meta data object with some desired attributes.
		vtkNew<vtkDICOMMetaData> meta; // !!! NOVO
		// !!! auto generator = vtkSmartPointer <vtkDICOMCTGenerator>::New(); // !!! VELHO
		// !!! vtkSmartPointer <vtkDICOMMetaData> meta = vtkSmartPointer <vtkDICOMMetaData>::New(); // !!! VELHO

		// VERIFICAR COMO LIDAR QUANDO AS TOMOGRAFIAS DE ENTRADA TEM ESCALAS DIFERENTES
		// !!! Quando le, troca os indices 0 e 1 !!!
		spacingDataO[0] = xScaleORG; // @@@
		spacingDataO[1] = yScaleORG; // @@@
		spacingDataO[2] = zScaleORG; // @@@

		cout << "spacingDataO: " << spacingDataO[0] << " , " << spacingDataO[1] << " , " << spacingDataO[2] << endl;

		int width = dimsO[0] - x0 + 1;
		int height = dimsO[1] - y0 + 1;
		int depth = dimsO[2] - z0 + 1;

		// Ajusta dimensoes do arquivo
		img->SetOrigin(0, 0, 0);
		img->SetDimensions(dimsO[0], dimsO[1], dimsO[2]); // @@@	
		img->SetSpacing(spacingDataO[0], spacingDataO[1], spacingDataO[2]); // @@@

		img->AllocateScalars(VTK_SHORT, 1);

		std::vector<uint16_t> raw_data(width*height*depth, 0);

		const double WIarr[4][4] = {
			{WI(0,0),WI(0,1),WI(0,2),WI(0,3)},
			{WI(1,0),WI(1,1),WI(1,2),WI(1,3)},
			{WI(2,0),WI(2,1),WI(2,2),WI(2,3)},
			{WI(3,0),WI(3,1),WI(3,2),WI(3,3)}
		};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Assembling tomography
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp parallel for schedule(guided)
		int dC = -1;
		for (int z = z0; z < dimsO[2]; z++)
		{
			int xlS[4], ylS[4], zlS[4];
			double p[4][4][4];

			double xiS, yiS, ziS;

			double wS, wT, cS, cT;

			for (int y = y0; y < dimsO[1]; y++)
			{
				for (int x = x0; x < dimsO[0]; x++)
				{
					// Inicializa pesos dos voxels: [S]ource e [T]arget
					wS = 1.0;
					wT = 1.0;

					// Inicializa valores efetivos dos voxels
					cS = 0.0;
					cT = 0.0;

					// Transformacao homogenea inversa para achar o ponto respectivo no volume original
					xiS = floor(WIarr[0][0] * x + WIarr[0][1] * y + WIarr[0][2] * z + WIarr[0][3]); // ###
					yiS = floor(WIarr[1][0] * x + WIarr[1][1] * y + WIarr[1][2] * z + WIarr[1][3]); // ###
					ziS = floor(WIarr[2][0] * x + WIarr[2][1] * y + WIarr[2][2] * z + WIarr[2][3]); // ###

					// Bordas nao sao confiaveis, sendo consideradas com um peso menor no blending
					if ((xiS < border) || (xiS >= dimsS[0] - border) || (yiS < border) || (yiS >= dimsS[1] - border) || (ziS < border) || (ziS > dimsS[2] - border))
						wS = 0.1;
					if ((x < border) || (x >= dimsT[0] - border) || (y < border) || (y >= dimsT[1] - border) || (z < border) || (z > dimsT[2] - border))
						wT = 0.1;

					// Verifica se as coordenadas sao validas (dentro dos limites da tomografia original)
					if ((xiS >= 0) && (xiS < dimsS[0]) && (yiS >= 0) && (yiS < dimsS[1]) && (ziS >= 0) && (ziS < dimsS[2])) 
						{

						//////////////////////////////////////////////////////////////////
						PERC_OLD = PERC;
						PERC = floor( 100.0 * ziS / (float)(dimsO[2]-1) ); // !!!
                		PERC = max(PERC_OLD, PERC + PERC_ACM); // !!!
						//////////////////////////////////////////////////////////////////

						cS = raw_dataS[ziS * heightS * widthS + yiS * widthS + xiS]; // @@@
						// cS = raw_dataS[int(ziS) * heightS * widthS + int(yiS) * widthS + int(xiS)]; // @@@
						// cS = tricubicInterpolate(p, xfrac, yfrac, zfrac); // @@@
					}
					else
						wS = 0.0;

					if ((x >= 0) && (x < dimsT[0]) && (y >= 0) && (y < dimsT[1]) && (z >= 0) && (z < dimsT[2]))
						{
						cT = raw_dataT[z * heightT * widthT + y * widthT + x];
						}
					else
						wT = 0.0;

					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					if (wT < 1e-10) 
						cT = 0;
					if (wS < 1e-10) 
						cS = 0;
					
					if (dC > 0)
						{
						if ( ((minzSo - dC) > (minzT + dC)) && ((minzSo + dC) < (maxzT-dC)) )
							{
							// s embaixo de t
							if (z > (maxzT - dC))
								wS = 1;
							if (z < (minzSo + dC))
								wS = 0;
							if (z <= (maxzT - dC) && z >= (minzSo + dC))
								{
								wS = ((z - (minzSo + dC)) / ((maxzT -dC) - (minzSo + dC))); // !!!
								}
							}
						else
							{
							// s em cima de t
							if (z < (minzT + dC))
								wS = 1;
							if (z > (maxzSo - dC))
								wS = 0;
							if (z >= minzT && z <= maxzSo)
								if (z >= (minzT + dC) && z <= (maxzSo -dC))
									{
									wS = 1 - ((z - (minzT + dC)) / ((maxzSo - dC) - (minzT + dC))); // !!!
									}
							}
						}

					wT = 1 - wS;

					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					// Blending
					double cO = (wS * cS + wT * cT) / (wS + wT); // !!!

					// Armazena voxel resultante:
					int xt, yt, zt;

					xt = x - x0;
					yt = y - y0;
					zt = z - z0;

					// Verificacao adicional de limites da matriz
					if ((xt < 0) || (xt >= dimsO[0]) || (yt < 0) || (yt >= dimsO[1]) || (zt < 0) || (zt >= dimsO[2]))
						continue;

					if (cO < 0) 
						cO = 0;

					raw_data[zt * height * width + yt * width + xt] = (uint16_t)std::lround(cO); // !!!

				}
			}

			//////////////////////////////////////////////////////////////////////////////////////////////
			// Deteccao de dC (slice inicial do cilindro em tomografia conica)
			//////////////////////////////////////////////////////////////////////////////////////////////
			double midm = floor(dimsO[0] / 2);
			double midn = floor(dimsO[1] / 2);
			double rr;
			for (int i = dimsO[0] - 1; i >= 0; i--)
			{
				if (raw_data[z * height * width + midn * width + i] > 0)
				{
					rr = abs((double)i - midm);
					break;
				}
			}

			/*
			// Desenha circunferencias com os raios dos slices
			for (double alpha = 0; alpha <= 2 * M_PI; alpha += 0.001)
			{
				double xx = floor(rr * cos(alpha)) + midm;
				double yy = floor(rr * sin(alpha)) + midn;
				if ((xx >= 0) && (xx < dimsO[0]) && (yy >= 0) && (yy < dimsO[2]))
					raw_data[z * height * width + yy * width + xx] = 4095;
			}
			*/

			if ( (dC == - 1) && (rr >= 300) )
			{
				dC = z;
				cout << endl << ">>>>>>>>> Profundidade z e raio rr: " << z << " | " << rr << endl;
			}

			// cout << endl << ">>>>>>>>> Profundidade z e raio rr: " << z << " | " << rr << endl;

		}

		cout << endl << ">>> Progress: " << PERC << endl << endl;
		PERC_ACM = 300;
		// PERC_ACM = PERC; // !!!

		for (int z = 0; z < dimsO[2]; z++)
		{
			PERC = floor(100.0 * z / (float)dimsO[2]);
			PERC = PERC + PERC_ACM;

			for (int y = 0; y < dimsO[1]; y++)
			{
				for (int x = 0; x < dimsO[0]; x++)
				{
					short int* pixel = static_cast<short int*>(img->GetScalarPointer(x, y, z));
					pixel[0] = (short int)((double)raw_data[z * height * width + y * width + x] + minden);
				}
			}
		}
		cout << "TRANSFORM ENDS" << endl;

		cout << endl << ">>> Progress: " << PERC << endl << endl;
 		PERC_ACM = PERC;

		cout << endl << "==========================================================================" << endl << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Grava arquivo raw_data
		stringstream filename;
		filename << file_out << "/raw_data_" << width << "x" << height << "_" << depth << ".raw";
		fstream file;
		file.open(filename.str(), ios::out | ios::binary);
		file.write((const char*)raw_data.data(), raw_data.size() * 2);
		file.close();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Escreve arquivo DICOM contendo a fusao
		vtkSmartPointer <vtkDICOMWriter> dicom_writer = vtkSmartPointer <vtkDICOMWriter>::New();

		// Configura arquivi DICOM
		dicom_writer->SetInputData(img);
		dicom_writer->SetMetaData(meta); // !!!
		dicom_writer->SetGenerator(generator); // !!!

		// dicom_writer->SetSeriesDescription("Sagittal Multi-planar Reformat");

		cout << "\n\nWriting DICOM file\n";
		// cout << "Dimensions: " << dimsO[0] << " , " << dimsO[1] << " , " << dimsO[2] << endl;

		// Set the output filename format as a printf-style string.
		dicom_writer->SetFilePattern("%s/%03d.dcm");

		// Set the directory to write the files into.
		dicom_writer->SetFilePrefix(file_out);

		// Write the file.
		dicom_writer->Write();

		cout << "DICOM Write Progress: 100%" << std::endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Le arquivo DICOM com a fusao
		cout << "\nReading DICOM FUSION" << std::endl;

		std::shared_ptr<geometry::PointCloud> fusion;

		fusion = readDicomToOpen3D(file_out, 255, 215, 0);
		fusion->EstimateNormals(open3d::geometry::KDTreeSearchParamHybrid(4, 30));

		cout << "DICOM Read Progress: 100%" << std::endl;

		// Exibe fusao
		// if (visualize)
		VisualizeComponent(*fusion);

		return 0;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Pre-processamento de pointcloud (Open3D) a partir de pointcloud (Open3D) 
	std::tuple<std::shared_ptr<geometry::PointCloud>, std::shared_ptr<registration::Feature>> PreprocessPointCloud(std::shared_ptr<geometry::PointCloud> pcd)
	{
		// Pre-processamento de pointcloud (Open3D): preprocess a partir da memoria

		// Faz downsample
		auto pcd_down = pcd->VoxelDownSample(decimation); //  2 1 //  Depois que aumentei para valor acima de 0.6, conseguiu eliminar vertices
		// !!! auto pcd_down = pcd->VoxelDownSample(decimation); //  2 1 //  Depois que aumentei para valor acima de 0.6, conseguiu eliminar vertices // @@@

		// Estima as normais depois do downsample
		pcd_down->EstimateNormals(open3d::geometry::KDTreeSearchParamHybrid(2 * decimation, 30)); //  4 2 //  Com o valor 1.2, conseguiu estimar normais (2 vezes maior do que o VoxelDownSample 0.6)

		//  Opcional
		pcd_down->NormalizeNormals();

		// Calcula vetor de features FPFH33 (vetor FPFH com 33 elementos)
		auto pcd_fpfh = registration::ComputeFPFHFeature(*pcd_down, open3d::geometry::KDTreeSearchParamHybrid(10 * decimation, 100)); // 20,100  Funcionou com 20,100 para downsample 2 e estimatenormals 4,30

		return std::make_tuple(pcd_down, pcd_fpfh);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void VisualizeComponent(const open3d::geometry::PointCloud& source)
	{
		std::shared_ptr<geometry::PointCloud> source_ptr(new geometry::PointCloud);
		*source_ptr = source;

		// Transformacao para visualizacao (rotacao):
		Eigen::Affine3d rot;
		rot = Eigen::AngleAxisd(M_PI / 2.0, Eigen::Vector3d(1.0, 0.0, 0.0));
		Eigen::Matrix4d vis4d;
		vis4d = rot.matrix();
		source_ptr->Transform(vis4d);

		visualization::DrawGeometries({ source_ptr }, "Component - Type q to proceed", 1600, 1200);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void VisualizeBeforeRegistration(const open3d::geometry::PointCloud& source,
		const open3d::geometry::PointCloud& target)
	{
		std::shared_ptr<geometry::PointCloud> source_ptr(new geometry::PointCloud);
		std::shared_ptr<geometry::PointCloud> target_ptr(new geometry::PointCloud);
		*source_ptr = source;
		*target_ptr = target;

		// Transformacao para visualizacao (rotacao e translacao):
		Eigen::Affine3d rot;
		rot = Eigen::AngleAxisd(M_PI / 2.0, Eigen::Vector3d(1.0, 0.0, 0.0));
		Eigen::Affine3d tra;
		tra = Eigen::Translation3d(Eigen::Vector3d(0.0, 0.0, 0.0));
		Eigen::Affine3d vis3d;
		vis3d = rot * tra;
		Eigen::Matrix4d vis4d;
		vis4d = vis3d.matrix();
		source_ptr->Transform(vis4d);
		tra = Eigen::Translation3d(Eigen::Vector3d(0.0, 0.0, 260.0)); // Translacao inicial da componente (apenas para efeito visual)		
		// tra = Eigen::Translation3d(Eigen::Vector3d(0.0, 0.0, 260.0)); // Translacao inicial da componente (apenas para efeito visual) // @@@		
		vis3d = rot * tra;
		vis4d = vis3d.matrix();
		target_ptr->Transform(vis4d);

		visualization::DrawGeometries({ source_ptr, target_ptr }, "Result without registration - Type q to proceed", 1600, 1200);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void VisualizeAfterRegistration(const open3d::geometry::PointCloud& source,
		const open3d::geometry::PointCloud& target,
		const Eigen::Matrix4d& Transformation)
	{
		std::shared_ptr<geometry::PointCloud> source_transformed_ptr(new geometry::PointCloud);
		std::shared_ptr<geometry::PointCloud> target_ptr(new geometry::PointCloud);
		*source_transformed_ptr = source;
		*target_ptr = target;

		// Transformacao para registration (stitching)
		source_transformed_ptr->Transform(Transformation);

		// Transformacao para visualizacao (rotacao pura):
		Eigen::Affine3d rot;
		rot = Eigen::AngleAxisd(M_PI / 2.0, Eigen::Vector3d(1.0, 0.0, 0.0));
		Eigen::Matrix4d vis4d;
		vis4d = rot.matrix();
		source_transformed_ptr->Transform(vis4d);
		target_ptr->Transform(vis4d);

		visualization::DrawGeometries({ source_transformed_ptr, target_ptr }, "Registration result - Type q to proceed", 1600, 1200);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Le arquivo DICOM para pointcloud (Open3D)
	// !!! LER DA MEMORIA !!!
	std::shared_ptr<geometry::PointCloud> readDicomToOpen3DMemory(std::shared_ptr<geometry::PointCloud> my_cloud_ptr)
	{
		auto cloud_ptr = std::make_shared<geometry::PointCloud>(); //  Open3D
		
		cloud_ptr = my_cloud_ptr; //  Open3D

		return cloud_ptr;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Grava arquivo DICOM para pointcloud (Open3D)
	// !!! GRAVAR NA MEMORIA !!!
	void writeDicomToOpen3DMemory(std::shared_ptr<geometry::PointCloud> my_cloud_ptr, std::shared_ptr<geometry::PointCloud> cloud_ptr)
	{		
		my_cloud_ptr = cloud_ptr; //  Open3D
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Le arquivo DICOM para pointcloud (Open3D)
	std::shared_ptr<geometry::PointCloud> readDicomToOpen3D(const char *dcmfolderpath, unsigned char r, unsigned char g, unsigned char b)
	{
		pcl::PointCloud<pcl::POINT_TYPE>::Ptr mx_cloud_ptr(new pcl::PointCloud<pcl::POINT_TYPE>); //  PCL

		auto cloud_ptr = std::make_shared<geometry::PointCloud>(); //  Open3D

		pcl::PointCloud<pcl::POINT_TYPE>::Ptr cloud(new pcl::PointCloud<pcl::POINT_TYPE>);

		// DICOM read
		vtkFileOutputWindow* outwin = vtkFileOutputWindow::New();
		outwin->SetFileName("vtklog.txt");
		outwin->SetInstance(outwin);

		vtkSmartPointer<ErrorObserver> errorObserver = vtkSmartPointer<ErrorObserver>::New();
		vtkSmartPointer<vtkDICOMImageReader> reader = vtkSmartPointer<vtkDICOMImageReader>::New();

		reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
		reader->AddObserver(vtkCommand::WarningEvent, errorObserver);

		sliceData = vtkSmartPointer<vtkImageData>::New();

		reader->SetDirectoryName(dcmfolderpath); // dcmfolderpath
		reader->Update();

		if (errorObserver->GetError())
		{
			std::cout << "Failed DICOM file access: " << std::endl;
			std::cout << errorObserver->GetErrorMessage();

			pcl::PointCloud<pcl::POINT_TYPE>::Ptr pt; //
			return NULL; // Nao compila em alguns ambientes
		}

		sliceData = reader->GetOutput();
		outwin->Delete();

		// Extracting pointcloud
		int numberOfDims = sliceData->GetDataDimension();

		dims = sliceData->GetDimensions();
		std::cout << "Cloud dimensions: ";
		int totalPoints = 1;
		for (int i = 0; i < numberOfDims; i++)
		{
			std::cout << dims[i] << " , ";
			totalPoints = totalPoints * dims[i];
		}
		std::cout << std::endl;

		std::cout << "Reading DICOM file" << std::endl;
		std::cout << "Number of DICOM points: " << totalPoints << std::endl;

		// Generate PCL point cloud
		double* dataRange = sliceData->GetScalarRange();
		double* spacingData = reader->GetDataSpacing();

		std::cout << "Data intensity bounds... min: " << dataRange[0] << ", max: " << dataRange[1] << std::endl;
		if (numberOfDims != 3)
		{
			std::cout << "Incorrect number of dimensions in DICOM file, generation failed..." << std::endl;

			pcl::PointCloud<pcl::POINT_TYPE>::Ptr pt;
			return NULL; // Nao compila em alguns ambientes
		}
		else
		{
			// !!! Scaling a partir do arquivo DICOM (correto)
			xScaleORG = spacingData[0]; // @@@
			yScaleORG = spacingData[1]; // @@@
			zScaleORG = spacingData[2]; // @@@

			cout << ">>> Original dimensions from DICOM file: " << xScaleORG << ", " << yScaleORG << ", " << zScaleORG << endl;

			//  !!! Ignorei scaling
			xScale = 1.0; //  !!! Ignorei scaling
			yScale = 1.0; //  !!! Ignorei scaling
			zScale = 1.0; //  !!! Ignorei scaling

			cout << ">>> Forced dimensions to DICOM file: " << xScale << ", " << yScale << ", " << zScale << endl;

			Scale = spacingData[0];

			std::cout << "x spacing: " << xScale << std::endl;
			std::cout << "y spacing: " << yScale << std::endl;
			std::cout << "z spacing: " << zScale << std::endl;

			for (int z = 0; z < dims[2]; z++)
			{
				for (int x = 0; x < dims[0]; x++)
				{
					for (int y = dims[1] - 1; y >= 0; y--) //  De traz para a frente
					{
						pcl::POINT_TYPE tempPt = pcl::POINT_TYPE();

						tempPt.x = x * xScale;
						tempPt.y = y * yScale;
						tempPt.z = z * zScale;

						// Vetor normal (Recalculado no programa principal)
						// Esses valores iniciais para fora da superficie sao importantes para evitar que
						// sejam estimadas normais para dentro do cranio, posteriormente na estimativa das normais
						// Esses valores precisam ser ajustados em caso de rotacao para visualizacao
						tempPt.normal_x = 0; // 0; 0 (com rotacao de 90 graus)
						tempPt.normal_y = 1; // 1; 0 (com rotacao de 90 graus)
						tempPt.normal_z = 0; // 0; 1 (com rotacao de 90 graus)

						// Cor RGB (precisa preencher novamente no programa principal)
						tempPt.r = 0;
						tempPt.g = 255;
						tempPt.b = 0;

						double tempIntensity = sliceData->GetScalarComponentAsDouble(x, y, z, 0);

						if (isinf(tempIntensity))
							tempIntensity = -1;

						//  O limiar para extrair o osso esta aqui (1376)
						if (tempIntensity >= thr_bone) //  1376
						{
							cloud->points.push_back(tempPt);
							break; //  Pego apenas o primeiro ponto (superficie/casca)
						}
					}
				}
				//std::cout << reader->GetTransferSyntaxUID() << std::endl;
				//std::cout << reader->get
			}

		}

		std::cout << "DICOM Read Progress: 100%" << std::endl;

		mx_cloud_ptr = cloud;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		size_t num_points = mx_cloud_ptr->size();
		printf("Size of pointcloud read from DICOM file (thresholded): %zu\n", num_points);

		// Alocacao de memoria
		cloud_ptr->points_.resize(num_points);
		cloud_ptr->normals_.resize(num_points);
		cloud_ptr->colors_.resize(num_points);

		// Conversao de pointcloud PCL para pointcloud Open3D
		for (size_t i = 0; i < num_points; i++)
		{
			cloud_ptr->points_[i](0) = mx_cloud_ptr->points[i].x;
			cloud_ptr->points_[i](1) = mx_cloud_ptr->points[i].y;
			cloud_ptr->points_[i](2) = mx_cloud_ptr->points[i].z;

			cloud_ptr->normals_[i](0) = mx_cloud_ptr->points[i].normal_x;
			cloud_ptr->normals_[i](1) = mx_cloud_ptr->points[i].normal_y;
			cloud_ptr->normals_[i](2) = mx_cloud_ptr->points[i].normal_z;

			cloud_ptr->colors_[i](0) = (float)r / (float)255;
			cloud_ptr->colors_[i](1) = (float)g / (float)255;
			cloud_ptr->colors_[i](2) = (float)b / (float)255;
		}

		utility::LogInfo("End of the test.\n");

		return cloud_ptr;
	}

}; // Fim da classe

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int stitch_main(int argc, char* argv[])
{
	// Progress bar
	// std::thread my_bar_text(progress_bar_text); // spawn new thread
	std::thread my_bar_graphic(progress_bar_graphic); // spawn new thread

	// Faz o stitching de duas tomografias DICOM, resultando em uma terceira tomografia DICOM
	if (argc != 4)
	{
		utility::LogInfo("Usage : Testvisualizer [path_to_first_DICOM_folder_input] [path_to_second_DICOM_folder_input] [path_to_third_DICOM_folder_output]");
		return 1;
	}

     // Cria objeto Stitcher
    Stitcher S;
        
    // Tomografias: FR OR MX MD
	S.stitch(argv[1], argv[2], argv[3]);
	// S.stitch("../dataset/0001_ANONIMO_ANONIMO_8X16_MX_2019-08-21", "../dataset/0001_ANONIMO_ANONIMO_8X16_MD_2019-08-21", "../dataset/fusion1");

	// synchronize threads:
	// my_bar_text.join();
	my_bar_graphic.join();

    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 