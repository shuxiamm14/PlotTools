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
	fstream *file;
	(*file).open(filename.c_str());
	(*file)<<"\\footnotesize\n";
	(*file)<<"\\caption{"<<caption<<"}\n";
	(*file)<<"\\centering\n";
	int ncolumn = columns.size();
	if(ncolumn > maxcolumn){
		int nlong = ncolumn%(ncolumn/maxcolumn+1);	   //10%3 = 1
		int count = 0;
		int averagelow = ncolumn/(ncolumn/maxcolumn+1); // 10/3 = 3

		vector<string> new_columns;
		for (int ivec = 0; ivec < count; ++ivec)
		{
			for (int i = 0; i < averagelow+1; ++i)
			{
				if(ivec == averagelow && ivec >= nlong) continue;
				new_columns.push_back(columns[count]);
				count ++;
			}
			writeContent(new_columns, &(*file));
			new_columns.clear();
		}
	}
	(*file)<<"\\label{tab:"<<label<<"}\n";
	(*file)<<"\\end{table}\n";
}

void LatexChart::writeContent(std::vector<std::string> new_columns, std::fstream* file){
	(*file)<<"\\begin{tabular}{|";
	for(auto new_column: new_columns) (*file)<<"c|";
	(*file)<<"} \\hline\n";
	//==============================column title=====================================
	for(auto new_column: new_columns) (*file)<<" & "<<new_column;
	(*file)<<"\\\\\\hline\n";
	//==============================table content=====================================
	for(auto row: rows){
		(*file)<<row;
		for(auto new_column: new_columns) (*file)<<" & "<<content[row][new_column].nominal<<"+/-"<<content[row][new_column].error;
		(*file)<<"\\\\\\hline\n";
	}
	(*file)<<"\\end{tabular}\n";
}