#include "makechart.h"
#include <string>
#include <fstream>
#include <algorithm>
using namespace std;

void LatexChart::set(std::string row, std::string column, float nominal, float error){
	observable obs(nominal,error);
	set(row,column,obs);
}

void LatexChart::set(std::string row, std::string column, double nominal, double error){
	observable obs(nominal,error);
	set(row,column,obs);
}

void LatexChart::set(std::string row, std::string column, observable obs){
	auto iter = find(rows.begin(),rows.end(),row);
	if (iter == rows.end()) rows.push_back(row);
	iter = find(columns.begin(),columns.end(),column);
	if (iter == columns.end()) columns.push_back(column);
	content[row][column] = obs;
}

void LatexChart::clear(){
	rows.clear();
	columns.clear();
	content.clear();
}

void LatexChart::reset(){
	for(auto row : content)
		for(auto column : row.second)
			column.second = 0;
}

void LatexChart::print(std::string filename){
	fstream file;
	file.open(filename.c_str());
	file<<"\\footnotesize\n";
	file<<"\\caption{"<<caption<<"}\n";
	file<<"\\centering\n";
	file<<"\\begin{tabular}{|";
	for(auto column: columns) file<<"c|";
	file<<"} \\hline\n";
	//==============================column title=====================================
	for(auto column: columns) file<<" & "<<column<<"\n";
	file<<"\\\\\\hline\n";
	//==============================table content=====================================
	for(auto row: rows){
		file<<row;
		for(auto column: columns) file<<" & "<<content[row][column].nominal<<"+/-"<<content[row][column].error<<"\n";
		file<<"\\\\\\hline\n";
	}
	//==============================end table=====================================
	file<<"\\end{tabular}\n";
	file<<"\\label{tab:"<<label<<"}\n";
	file<<"\\end{table}\n";
}
