#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

class CooTrans
{
public:
	double iPI = 0.0174532925199433; //3.1415926535898/180.0
	double PI = 3.1415926535898;
	// 54年北京坐标系参数
	// double a = 6378245.0;   // 长轴
	// double f = 1.0 / 298.3; // 扁率   (a-b)/a

	// 80年西安坐标系参数
	// double a=6378140.0;
	// double f=1/298.257;

	// WGS84坐标系参数
	double a = 6378137.0;
	double f = 1 / 298.257223563;

	double ZoneWide = 6.0;     // 带宽
	double e2 = 2 * f - f * f; // e为第一偏心率，可以算可以直接提供，e2 = e * e

	//正算 输入为经纬度，并非是分秒格式，输出为高斯坐标
	std::vector<double> LLtoGaussXY(double longitude, double latitude)
	{
		int ProjNo = (int)(longitude / ZoneWide);
		double longitude0 = ProjNo * ZoneWide + ZoneWide / 2;
		longitude0 = longitude0 * iPI;

		// 经度转换为弧度
		longitude = longitude * iPI;
		// 纬度转换为弧度
		latitude = latitude * iPI;

		double ee = e2 * (1.0 - e2);
		double N = a / sqrt(1.0 - e2 * sin(latitude) * sin(latitude)); // 该点的卯酉圈曲率半径
		double T = tan(latitude) * tan(latitude);
		double C = ee * cos(latitude) * cos(latitude);
		double A = (longitude - longitude0) * cos(latitude);
		double M = a * ((1 - e2 / 4 - 3 * e2 * e2 / 64 - 5 * e2 * e2 * e2 / 256) * latitude - (3 * e2 / 8 + 3 * e2 * e2 / 32 + 45 * e2 * e2 * e2 / 1024) * sin(2 * latitude) + (15 * e2 * e2 / 256 + 45 * e2 * e2 * e2 / 1024) * sin(4 * latitude) - (35 * e2 * e2 * e2 / 3072) * sin(6 * latitude));
		double xval = N * (A + (1 - T + C) * A * A * A / 6 + (5 - 18 * T + T * T + 72 * C - 58 * ee) * A * A * A * A * A / 120);
		double yval = M + N * tan(latitude) * (A * A / 2 + (5 - T + 9 * C + 4 * C * C) * A * A * A * A / 24 + (61 - 58 * T + T * T + 600 * C - 330 * ee) * A * A * A * A * A * A / 720);
		double X0 = 500000;
		double Y0 = 0;
		xval = xval + X0;
		yval = yval + Y0;
		return { xval, yval };
	}
	//反算 输入为高斯坐标，输出为经纬度
	std::vector<double> GaussXYtoLL(double X, double Y)
	{
		int ProjNo = (int)(X / 1000000); // 查找带号
		double longitude0 = (ProjNo - 1) * ZoneWide + ZoneWide / 2;
		longitude0 = longitude0 * iPI; // 中央经线

		double X0 = ProjNo * 1000000 + 500000;
		double Y0 = 0;
		double xval = X - X0;
		double yval = Y - Y0; // 带内大地坐标

		double e1 = (1.0 - sqrt(1 - e2)) / (1.0 + sqrt(1 - e2));
		double ee = e2 / (1 - e2);
		double M = yval;
		double u = M / (a * (1 - e2 / 4 - 3 * e2 * e2 / 64 - 5 * e2 * e2 * e2 / 256));
		double fai = u + (3 * e1 / 2 - 27 * e1 * e1 * e1 / 32) * sin(2 * u) + (21 * e1 * e1 / 16 - 55 * e1 * e1 * e1 * e1 / 32) * sin(4 * u) + (151 * e1 * e1 * e1 / 96) * sin(6 * u) + (1097 * e1 * e1 * e1 * e1 / 512) * sin(8 * u);
		double C = ee * cos(fai) * cos(fai);
		double T = tan(fai) * tan(fai);
		double N = a / sqrt(1.0 - e2 * sin(fai) * sin(fai)); // 该点的卯酉圈曲率半径
		double R = a * (1 - e2) / sqrt((1 - e2 * sin(fai) * sin(fai)) * (1 - e2 * sin(fai) * sin(fai)) * (1 - e2 * sin(fai) * sin(fai)));
		double D = xval / N;
		// 计算经度(Longitude) 纬度(Latitude)
		double longitude = longitude0 + (D - (1 + 2 * T + C) * D * D * D / 6 + (5 - 2 * C + 28 * T - 3 * C * C + 8 * ee + 24 * T * T) * D * D * D * D * D / 120) / cos(fai);
		double latitude = fai - (N * tan(fai) / R) * (D * D / 2 - (5 + 3 * T + 10 * C - 4 * C * C - 9 * ee) * D * D * D * D / 24 + (61 + 90 * T + 298 * C + 45 * T * T - 256 * ee - 3 * C * C) * D * D * D * D * D * D / 720);
		// 转换为度 DD
		return { longitude / iPI, latitude / iPI };
	}
};

CooTrans ct;

double getDegree(double degreeMinute) {
	double res = 0;
	double integerPos = floor(degreeMinute);
	double littlePos = degreeMinute - integerPos;
	res = littlePos * 100 / 60 + integerPos;
	return res;
}

// 118.69580148,32.133512018999987,64.930999999358178
double lon_orgin = 118.69580148;
double lat_origin = 32.133512018999987;
double alt_origin = 64.930999999358178; // 坐标增加的数值

int main (int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input.pcd> <output.pcd>" << std::endl;
        return -1;
    }

    // 创建一个 PointCloud<pcl::PointXYZRGB> 对象，用于存储输入和输出的点云
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGB>);

    // 读取点云文件
    if (pcl::io::loadPCDFile<pcl::PointXYZRGB> (argv[1], *cloud) == -1) {
        PCL_ERROR ("Couldn't read file %s\n", argv[1]);
        return -1;
    }

    std::cout << "Loaded " << cloud->width * cloud->height
              << " data points from " << argv[1] << std::endl;

    // 处理每个点
    for (auto& point : *cloud) {

        point.x += lon_orgin;  // 增加 x 坐标
        point.y += lat_origin;  // 增加 y 坐标
        point.z += alt_origin;  // 增加 z 坐标

        double lon = point.x , lat= point.y , alt=point.z;
        lon = getDegree(lon);
        lat = getDegree(lat);

        double x, y, z;
        std::vector<double> coordinate = ct.LLtoGaussXY(lon, lat);
        x = coordinate.at(1);
        y = coordinate.at(0);
        z = alt;

		point.x = x;
		point.y = y;
		point.z = z;

    }

    // 保存点云到新的 PCD 文件
    pcl::io::savePCDFileASCII (argv[2], *cloud);
    std::cout << "Saved " << cloud->width * cloud->height
              << " data points to " << argv[2] << std::endl;

    return 0;
}