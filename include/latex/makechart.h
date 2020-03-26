#include "observable.h"
#include "iostream"
#include <fstream>
#include <map>
class LatexChart
{
public:
	LatexChart(){ maxcolumn = 4; }
	~LatexChart(){};
	LatexChart(std::string chartlabel) : label(chartlabel) { maxcolumn = 4; }
	int maxcolumn;
	std::string label;
	std::string caption;
	std::vector<std::string> rows;
	std::vector<std::string> columns;
	std::map<std::string, std::map<std::string, observable>> content;
	void set(std::string row, std::string column, float nominal = 0, float error = 0, float errordown = 0);
	void set(std::string row, std::string column, double nominal = 0, double error = 0, double errordown = 0);
	void set(std::string row, std::string column, observable obs);
	void clear();
	void reset();
	void print(std::string filename);
	void writeContent(std::vector<std::string> new_columns, std::ofstream* file);
};
