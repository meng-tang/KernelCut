#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <graph.h>
#include <EasyBMP.h>
#include "basicutil.h"

#include "ncutknnbinary.h"
#include "ncutknnmulti.h"

enum Method {NCUTKNNBINARY, NCUTKNNMULTI};
Method method;
enum UserInput {BOX, SEEDS, FROMIMG};
UserInput userinput = BOX;
enum ErrorMeasure {ERRORRATE, FMEASURE, JACCARD, NOMEASURE};
ErrorMeasure errormeasure = ERRORRATE;

int main(int argc, char * argv[])
{
    srand( (unsigned)time( NULL ) );
    double totaltime = 0; // timing
    const char * UsageStr = "Usage: main -d DBdirectory -i imagename -n numberofSegmengts [-o outputdirectory] [-s w_smoothness] [-e errormeasure] [-u userinput] [-h on or off (hardconstraints)] [-k knn_k knnfiledirectory]\n";
    bool hardconstraintsflag = true;
    char * dbdir = NULL, * imgname = NULL, * methodstr = NULL, * outputdir = NULL;
    char * knnfiledir = NULL;
    char * initlabelingimgpath = NULL;
    Table2D<int> knntable;
    double w_smooth = 0;
    int KNN_K, numSegments;
    if(argc == 1){ 
        printf("%s",UsageStr);
        exit(-1);
    }
    for(int i=0;i<argc;i++){
        printf("%s ",argv[i]);
    }
    printf("\n");
    
    int opt;
    while((opt = getopt(argc, argv, "dineuowshbxkpg")) != -1){
        switch(opt){
            case 'd':
                dbdir = argv[optind];
                break;
            case 'i':
                imgname = argv[optind];
                break;
            case 'n':
                numSegments = atoi(argv[optind]);
                if(2 == numSegments){
                    methodstr = "NCUTKNNBINARY";
                    method = NCUTKNNBINARY;
                }
                else if(2 < numSegments){
                    methodstr = "NCUTKNNMULTI";
                    method = NCUTKNNMULTI;
                }
                else{
                    printf("number of segments not valid!\n");
                    exit(-1);
                }
                break;
            case 'e':
                if(0 == strcmp(argv[optind],"errorrate"))
                    errormeasure = ERRORRATE;
                else if(0 == strcmp(argv[optind],"fmeasure"))
                    errormeasure = FMEASURE;
                else if(0 == strcmp(argv[optind],"jaccard"))
                    errormeasure = JACCARD;
                else if(0 == strcmp(argv[optind],"nomeasure"))
                    errormeasure = NOMEASURE;
                else{
                    printf("measure not valid!\n");
                    exit(-1);
                }
                break;
            case 'u':
                if(0 == strcmp(argv[optind],"box"))
                    userinput = BOX;
                else if(0 == strcmp(argv[optind],"seeds"))
                    userinput = SEEDS;
                else if(0 == strcmp(argv[optind],"fromimage")){
                    userinput = FROMIMG;
                    initlabelingimgpath = argv[optind+1];
                }
                else{
                    printf("User input not valid!\n");
                    exit(-1);
                }
                break;
            case 'o':
                outputdir = argv[optind];
                break;
            case 's':
                w_smooth = atof(argv[optind]);
                break;
            case 'h':
                if(0 == strcmp(argv[optind],"on"))
                    hardconstraintsflag = true;
                else if(0 == strcmp(argv[optind],"off")){
                    hardconstraintsflag = false;
                }
                break;
            case 'k':
                KNN_K = atoi(argv[optind]);
                knnfiledir = argv[optind+1];
                break;
            default: /* '?' */
                fprintf(stderr, "%s", UsageStr);
                break;
        }
    }
    if((dbdir != NULL) && (imgname != NULL) && (methodstr != NULL))
        printf("database %s\nimage %s\nmethod %s\n", dbdir, imgname, methodstr);
    else
        printf("%s",UsageStr);


    printf("image : %s\n",imgname);
        
    // Read the RGB image
    Image image = Image((dbdir+string("/images/")+string(imgname) + string(".bmp")).c_str(),imgname,16,8);
        
    int imgw = image.img_w;
    int imgh = image.img_h;
        
            
    // Initial labeling
    Table2D<Label> initlabeling(imgw,imgh,NONE);
    Table2D<Label> hardconstraints(imgw,imgh,NONE);
    if(BOX == userinput){
        initlabeling = getinitlabeling(loadImage<RGB>((dbdir+string("/boxes/")+string(imgname) + string(".bmp")).c_str()),0);
        for(int i=0;i<imgw;i++){
	    for(int j=0;j<imgh;j++){
	        if(initlabeling[i][j]==BKG) hardconstraints[i][j] = BKG;
		else hardconstraints[i][j] = NONE;
            }
        }
        //image.addboxsmooth(getROI(hardconstraints,NONE));
    }
    else if((SEEDS == userinput) && (method != NCUTKNNMULTI)){
        initlabeling = getinitlabelingFB(loadImage<RGB>((dbdir+string("/seeds/")+string(imgname) + string(".bmp")).c_str()), red, blue);
        hardconstraints = initlabeling;
    }
    else if(FROMIMG == userinput){
        Table2D<RGB> initlabelingimg = loadImage<RGB>(initlabelingimgpath);
        for(int i=0;i<imgw;i++){
            for(int j=0;j<imgh;j++){
	        if(initlabelingimg[i][j]==white)
		    initlabeling[i][j] = BKG;
	        else
	           initlabeling[i][j] = OBJ;
	    }
	}
        hardconstraints = initlabeling;
    }
        
    Table2D<int> initlabeling_multi(imgw,imgh,0); // for multilabel
    RGB colors[6] = {white,red,green,blue,black,navy};
    if(method == NCUTKNNMULTI){
        if(userinput == SEEDS)
            initlabeling_multi = getinitlabelingMULTI(loadImage<RGB>((dbdir+string("/seedsmulti/")+string(imgname) + string(".bmp")).c_str()), colors, numSegments);
        else if(userinput == FROMIMG)
            initlabeling_multi = getinitlabelingMULTI(loadImage<RGB>(initlabelingimgpath), colors, numSegments);
        hardconstraints = initlabeling_multi;
    }
        
    if(hardconstraintsflag==false) hardconstraints.reset(NONE);
        
    // read knn graph
    if(method == NCUTKNNBINARY || method == NCUTKNNMULTI){
        char knnfile[100] = {0};
        strcat(knnfile,knnfiledir);
        strcat(knnfile,"/");
        strcat(knnfile,imgname);
        strcat(knnfile,".bin");
        //printf("knn file path:%s\n",knnfile);
        //read knn table
        Table2D<int> temp_knntable; 
        readbinfile(temp_knntable,knnfile, KNN_K/*8*/,imgw * imgh);
        knntable = Table2D<int>(imgw*imgh,KNN_K); // index from zero, KNN_K rows
        for(int i=0;i<KNN_K;i++){
	    for(int j=0;j<imgw*imgh;j++){
	        knntable[j][i] = temp_knntable[i/*8*/][j]-1; // make index start from zero
	    }
        }
        temp_knntable.resize(1,1);
    }
            
    clock_t start = clock();
    Table2D<Label> solution;
    Table2D<int> solution_multi;
	   
    if(method == NCUTKNNBINARY){
        solution = ncutknnbinarysegmentation(image, knntable, w_smooth, initlabeling, hardconstraints);
        knntable.resize(1,1);
    }
    else if(method == NCUTKNNMULTI){
        solution_multi = ncutknnmultisegmentation(image, knntable, w_smooth, initlabeling_multi, numSegments);
        knntable.resize(1,1);
    }
    clock_t finish = clock();
    totaltime = (double)(finish-start)/CLOCKS_PER_SEC;
        
    // save output
    if((outputdir!=NULL) &&( method != NCUTKNNMULTI))
        //savebinarylabeling(image.img, solution,(outputdir+string("/")+string(imgname) +"_s"+pch::to_string(w_smooth)+string(".bmp")).c_str());
	savebinarylabeling(image.img, solution,(outputdir+string("/")+string(imgname) +string("_")+string(methodstr)+"_s"+pch::to_string(w_smooth)+string(".bmp")).c_str());
    if((outputdir!=NULL) &&( method == NCUTKNNMULTI))
	savemultilabeling(solution_multi,(outputdir+string("/")+string(imgname) +string("_")+string(methodstr)+"_s"+pch::to_string(w_smooth)+string(".bmp")).c_str(), colors,image.img);

    // measures
    if(NOMEASURE != errormeasure)
    {
        // ground truth
        Table2D<Label> gt = getinitlabeling(loadImage<RGB>((dbdir+string("/groundtruth/")+string(imgname) + string(".bmp")).c_str()),255,0);
        double measure;
        if(ERRORRATE == errormeasure){
	    if( userinput == BOX && hardconstraintsflag )
	        measure = geterrorcount(solution,gt)/ (double) countintable(hardconstraints,NONE);
	    else
	        measure = geterrorcount(solution,gt)/ (double)(imgw*imgh);
        }
        else if(JACCARD == errormeasure)
	    measure = jaccard(solution, gt, OBJ);
        else if(FMEASURE == errormeasure)
	    measure = fmeasure(solution, gt, OBJ);
        printf("measure %.3f\n",measure);
    }

    cout<<"time for segmentation "<<totaltime<<" seconds!"<<endl;
    return -1;
}
