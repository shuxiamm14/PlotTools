#include "makechart.h"
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>
using namespace std;

void LatexChart::set(std::string row, std::string column, float nominal, float errorup, float errordown){
	observable obs(nominal,errorup,errordown);
	set(row,column,obs);
}

void LatexChart::set(std::string row, std::string column, double nominal, double errorup, double errordown){
	observable obs(nominal,errorup,errordown);
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
	ofstream *file = new ofstream();
	(*file).open(filename+".tex");
	(*file)<<"\\begin{table}\n";
	(*file)<<"\\footnotesize\n";
	(*file)<<"\\caption{"<<caption<<"}\n";
	(*file)<<"\\centering\n";
	int ncolumn = columns.size();
	if(ncolumn > maxcolumn){
		int nvec = ncolumn/maxcolumn;
		if(ncolumn%maxcolumn) nvec+=1;
		int nlong = ncolumn%nvec;
		int averagelow = ncolumn/nvec;
		if(!ncolumn%maxcolumn) averagelow-=1;
		vector<string> new_columns;
		int count = 0;
		for (int ivec = 0; ivec < nvec; ++ivec)
		{
			for (int i = 0; i < averagelow+1; ++i)
			{
				if(i == averagelow && ivec >= nlong) continue;
				new_columns.push_back(columns[count]);
				count ++;
			}
			writeContent(new_columns, file);
			new_columns.clear();
		}
	}else{
		writeContent(columns, file);
	}
	(*file)<<"\\label{tab:"<<label<<"}\n";
	(*file)<<"\\end{table}\n";
	file->close();
}

void LatexChart::writeContent(std::vector<std::string> new_columns, std::ofstream* file){
	(*file)<<"\\begin{tabular}{|";
	for(auto new_column: new_columns) (*file)<<"c|";
	(*file)<<"c|} \\hline\n";
	//==============================column title=====================================
	for(auto new_column: new_columns) (*file)<<" & "<<new_column;
	(*file)<<"\\\\\\hline\n";
	//==============================table content=====================================
	(*file)<<fixed<<setprecision(2);
	for(auto row: rows){
		(*file)<<row;
		for(auto new_column: new_columns) {
			(*file)<<" & ";
			if(content[row].find(new_column) == content[row].end())
				(*file)<<" /";
			else{
				(*file)<<"$"<<content[row][new_column].nominal;
				if(content[row][new_column].error) {
					if(content[row][new_column].error == content[row][new_column].errordown)
						(*file)<<"\\pm"<<content[row][new_column].error;
					else
						(*file)<<"^{+"<<content[row][new_column].error<<"}_{-"<<content[row][new_column].errordown<<"}";
				}
				(*file)<<" $";
			}
		}
		(*file)<<"\\\\\\hline\n";
	}
	(*file)<<"\\end{tabular}\n";
}

void LatexChart::add(LatexChart *target){
	for(auto row: rows){
		for(auto column: columns){
			content[row][column] += target->content[row][column];
		}
	}
}

LatexChart* LatexChart::clone(){
	LatexChart *chart = new LatexChart(label);
	chart->maxcolumn = maxcolumn;
	chart->label = label;
	chart->caption = caption;
	chart->rows = rows;
	chart->columns = columns;
	chart->content = content;
	return chart;
};