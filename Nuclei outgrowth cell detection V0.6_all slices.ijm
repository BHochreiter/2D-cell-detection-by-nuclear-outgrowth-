// this macro generates a cell detection map using an outgrowth algorithm based on a nuclei-image
// written by Bernhard Hochreiter, 2021

	batch=false; //true or false : processes an entire folder automatically
	name=".tif"; //selects for specific names when batch processing

	render=true; //true or false : intermediate steps are displayed - makes process much slower
	
	closeImages=false; //only relevant if batch=false

//nuclei detection
	nuclei_channel_number=1;
	nuclei_min_size=20;
	thresh="Default";	//manual, Default,...

//outgrowth iterations:
	iter=5;

//cell detection map
	create_cell_map=true;
	color_mode="multi"; //single or multi
	saveImage=false;


//analyze cells
	analyze=true;


//do not change anything after this line
//############################################################

print("##########################");
print(currentTime()+" MACRO START");

if(batch==true){
	setBatchMode(true);
	source=getDirectory("Choose a Directory");
	fileList=getFileList(source);
	evalList=newArray(0);

	for (i = 0; i < fileList.length; i++) {
		if(indexOf(fileList[i], "cellmap") >= 0) {}else{
			if(indexOf(fileList[i], name) >= 0) {
				evalList=append(evalList,fileList[i]);
			}
		}
	}

	print(currentTime()+" BATCH mode ON -"+evalList.length+" files found");

	for (i = 0; i < evalList.length; i++) {
		open(source+evalList[i]);
		splitAndAnalyze();		
		run("Close All");
		x=round(((i+1)/evalList.length)*100);
		print(currentTime()+""+(i+1)+" of "+evalList.length+" images processed ("+x+"%)");
		print("-------------------------------");
	}
}
else{
	if(render==false){
	setBatchMode(true);
	}
	splitAndAnalyze();
	setSlice(nuclei_channel_number);	
	run("Select None");
	if(closeImages==true){
		run("Close All");
	}
}	
print(currentTime()+" MACRO FINISH");
print("##########################");

//###########################################################
function splitAndAnalyze() {
	getDimensions(width, height, channels, slices, frames);
	noZ=slices;

	for (p = 0; p < noZ; p++) {
		run("Duplicate...", "duplicate slices="+p+1);
		title=getTitle();
		outgrowth();
		close(title);
		close("cell_map");
	}

}
//###########################################################

function outgrowth() {
	//preamble
	run("Set Measurements...", "  redirect=None decimal=3");
	if (isOpen("Results")) {}else{run("Measure");run("Clear Results");}
	title=getTitle();
	imageDir=getDirectory("image");
	titleClean=replace(title, ".nd2", "");
	titleClean=replace(titleClean, ".tif", "");
	print(currentTime()+" processing "+title);
	setBackgroundColor(0, 0, 0);
	
	resultstart=nResults;
	
	if(roiManager("count")>0){
		roiManager("Deselect");
		roiManager("Delete");	
	}
	
	//generate nuclei map
	setSlice(nuclei_channel_number);
	run("Select None");
	run("Duplicate...", "title=[cellmap]");
	run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
	
	if(thresh=="manual"){
		run("Threshold...");
		waitForUser("please adjust threshold");
	}
	else{
	setAutoThreshold(thresh+" dark");
	}
	run("Convert to Mask");
	run("Watershed");
	
	
	selectWindow("cellmap");
	run("Analyze Particles...", "size="+nuclei_min_size+"-infinity add");
	
	nnuc=roiManager("count");
	
	getDimensions(width, height, channels, slices, frames);
	newImage("nucleimap", "8-bit black", width, height,nnuc);
	
	//generate nuclei map stack
	
	for (i = 0; i < nnuc; i++) {
		selectWindow("nucleimap");
		setSlice(i+1);
		roiManager("Select", i);
		setForegroundColor(255, 255, 255);
		run("Fill", "slice");
	}
	close("cellmap");
	
	//outgrowth iterations
	
	for (i = 0; i < iter; i++) {
		selectWindow("nucleimap");
		run("Select None");
		run("Dilate", "stack");
		run("Z Project...", "projection=[Sum Slices]");
		rename("doubleevents");
		setThreshold(1.0000, 256.0000);
		run("Convert to Mask");
		imageCalculator("Multiply stack", "nucleimap","doubleevents");
		close("doubleevents");
		prog=round((i+1)/iter*100);
		print(currentTime()+" Iteration "+i+1+" of "+iter+" "+loadingBar());
	}
	
	//add full cells to ROI manager

	print(currentTime()+" processing ROIs...");
	
	for (i = 0; i < nnuc; i++) {
		setSlice(i+1);
		setThreshold(127, 255);
		run("Create Selection");
		roiManager("Add");
		run("Select None");
	}
		
	//generate final images
	
	close("nucleimap");
	if(batch==false){
		setBatchMode(false);
	}
	if(create_cell_map==true){
		print(currentTime()+" generating pretty images...");
		getDimensions(width, height, channels, slices, frames);
		newImage("cell_map", "RGB black", width, height, 1);
		
		run("RGB Color");
		for (i = 0; i < nnuc; i++) {
			roiManager("Select",i+nnuc);
			if(color_mode=="single"){
				setForegroundColor(255, 255, 0); //cytosol colour	
			}
			if(color_mode=="multi"){
				a=round(random*200+55);
				b=round(random*200+55);
				c=round(random*200+55);
				setForegroundColor(round(a/2), round(b/2), round(c/2));
			}
		
			run("Fill", "slice");	
			roiManager("Select",i);
			if(color_mode=="single"){
				setForegroundColor(255, 0, 0); //nuclei colour	
			}
			if(color_mode=="multi"){
				setForegroundColor(a, b, c);
			}	
			run("Fill", "slice");
			roiManager("Select",i+nnuc);
			setForegroundColor(0, 0, 0); //border colour
			run("Draw", "slice");
		}
		if(saveImage==true){
			saveAs("TIFF", imageDir+titleClean+"_cellmap");
		}
	}
	
	
	if(analyze==true){
		print(currentTime()+" analyzing intensities...");
		selectWindow(title);
		run("Select None");
		roiManager("Show None");
	
	//create results Array
		
		getDimensions(width, height, channels, slices, frames);
		n=channels*slices*frames;
		sliceArray=newArray(n);
		for (i = 0; i < n; i++) {
			sliceArray[i]=i+1;
		}
	
	
		for (i = 0; i < nnuc; i++) {
			// measure whole cell	
			roiManager("Select", i+nnuc);
			getStatistics(area, mean, min, max);
			setResult("image", resultstart+i, title);
			setResult("slice", resultstart+i, p+1);
			setResult("ROI", resultstart+i, i+1);
			setResult("cellarea", resultstart+i, area);
			
			for (j = 0; j < n; j++) {
				setSlice(j+1);
				getStatistics(area, mean, min, max);
				setResult("Channel_"+sliceArray[j]+"_wholeCellMean", resultstart+i, mean);			
			}
			
			//measure nuclei
			roiManager("Select", i);
			getStatistics(area, mean, min, max);
			setResult("nucleiarea", resultstart+i, area);
					
			for (j = 0; j < n; j++) {
				setSlice(j+1);
				getStatistics(area, mean, min, max);
				setResult("Channel_"+sliceArray[j]+"_nucleiMean", resultstart+i, mean);			
			}
					
			//measure cytosol
			roiManager("Select", newArray(i,i+nnuc));
			roiManager("XOR");
			getStatistics(area, mean, min, max);
			setResult("cytosolarea", resultstart+i, area);
					
			for (j = 0; j < n; j++) {
				setSlice(j+1);
				getStatistics(area, mean, min, max);
				setResult("Channel_"+sliceArray[j]+"_cytosolMean", resultstart+i, mean);			
			}
		}
		print(currentTime()+" analyzed "+nnuc+" objects");
	}
}
//##################################################

function append(arr, value) {
	arr2 = newArray(arr.length+1);
	for (i=0; i<arr.length; i++)
		arr2[i] = arr[i];
	arr2[arr.length] = value;
	return arr2;
}

//#######################################################################################

function currentTime() {
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	if(hour<10){hour="0"+hour;}
	if(minute<10){minute="0"+minute;}
	if(second<10){second="0"+second;}
	time=""+hour+":"+minute+":"+second+" - ";
	return time;
}

//#########################################################

function loadingBar(){
	loaded=round(prog/5);
	nonloaded=20-loaded;
	bar="|";
	
	for (i = 0; i < loaded; i++) {
		bar=bar+"o";
	}
	for (i = 0; i < nonloaded; i++) {
		bar=bar+"#";
	}
	bar=bar+"|";
    return bar;
}
