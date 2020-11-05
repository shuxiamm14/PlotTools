#include "atlasstyle/AtlasCanvas.h"
#include "atlasstyle/AtlasLabels.h"
#include "atlasstyle/AtlasStyle.h"

atlasCanvas::atlasCanvas(const char* name, const char* title):TCanvas(name, title, 600, 600){
	SetAtlasStyle();
	cd();
	ATLASLabel(0.2,0.860,"Internal",kBlack, title);
}