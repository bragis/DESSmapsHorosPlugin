//
//  DESSmapsWindowController.m
//  DESSmaps
//
//  Created by Bragi Sveinsson on 1/16/16.
//
//

#include "DESSmapsWindowController.h"

@implementation DESSmapsWindowController


// We are using this initialization method instead of the more standard init() method.
- (id) initWithFilter:(DESSmapsFilter *)filter {
    self = [super initWithWindowNibName:@"DESSmapsWindow"];
    baseFilter = filter;
    
    
    // pixList is an array containing the slices in the exam. The viewerController can be thought of as the "Exam", I think (see http://www.osirix-viewer.com/OsiriXDevKeynote.html).
    // The syntax is [obj method:argument], where in C++ we would write obj->method(argument). So here we are first getting a viewerController object vc, then getting its pixList by calling vc->pixList(). We could have given an argument here to select series no. N in the 4D viewer by selecting the N-th pixList: vc->pixList(N). Giving no argument is equal to saying N=0.
    NSMutableArray *pixList = [[baseFilter viewerController] pixList];
    // DCMPix points to a slice (see http://www.osirix-viewer.com/OsiriXDevKeynote.html). We are selecting slice no. 0 (the first slice).
    DCMPix *curPix = [pixList objectAtIndex:0];
    
    // Get image size - assume all slices the same size.
    zsize = [pixList count];        // Number of pixels in z-dimension (i.e. slices)
    xsize = [curPix pwidth];        // Number of pixels in x-dimension
    ysize = [curPix pheight];       // Number of pixels in y-dimension
    ntpts = [[baseFilter viewerController] maxMovieIndex];      // Number of time ppoints (i.e. number of series in the 4D viewer)
    nxyzpts = xsize*ysize*zsize;    // Total number of points
    NSLog(@"The exam has dimensions (x,y,z,t) = (%li,%li,%li,%li).",xsize,ysize,zsize,ntpts);
    
    // Start by setting these all to false. One of them should be true though.
    twoSeparateDESSscans = false;
    twoCombinedDESSscans = false;
    singleDESSscan = false;
    if (ntpts == 4) {
        twoSeparateDESSscans = true;
    } else {
        twoSeparateDESSscans = false;
    }
    
    // Get the scan parameters (TR, TE, etc) from the Dicom header
    [self getFieldsFromDicom];
    
    return self;
}





// windowDidLoad is called when the window is fully loaded. Some explanations:
// init (or a substitute) is the first method that gets called. Initializes self.
// awakeFromNib is called after init. When the Nib is loaded, all objects are allocated and intialized, and then their outlets and actions are hooked up, then the Nib loader sends awakeFromNib to every object in the Nib. When the window is fully loaded, the windowDidLoad method is called. This why outlets can't be accessed in the initializer.
-(void) windowDidLoad
// Called after the GUI loads. Now we can start putting numbers into the text fields.
{
    
    //Simple message to show we made it.
    NSLog(@"Starting windowDidLoad function.");
    //NSRunInformationalAlertPanel(@"Info:",@"Image Calc Controller Loaded.", @"OK", nil, nil);
    
    //[self updateFieldsInGUI];
    //[threshholdField setStringValue:[NSString stringWithFormat:@"%g",threshhold]];
    
    // Testing the getFieldsFromDicom method.
    //[self getFieldsFromDicom];
}




// Get the various fields from the Dicom header.
- (void) getFieldsFromDicom {
    NSLog(@"Starting getFieldsFromDicom function.");
    float acq1flipAngle, acq2flipAngle, oneTouchAcq1flipAngle, acq1tr, acq1te, acq2tr, acq2te;
    ViewerController *viewer = [baseFilter viewerController];               // Select the 4D viewer containing the DESS scans.
    NSMutableArray *pixListAcq1 = [viewer pixList];      // Select the first series from the 4D viewer.
    // Select the first slice of the first series. Contained as array of DCMPix objects (pixels).
    DCMPix *pixAcq1 = [pixListAcq1 objectAtIndex:0];
    
    // Read another slice for testing
    printf("SliceInterval: %g\n", [pixAcq1 sliceInterval]);
    printf("SliceLocation: %g\n", [pixAcq1 sliceLocation]);
    printf("SliceThickness: %g\n", [pixAcq1 sliceThickness]);
    printf("SpacingBetweenSlices: %g\n", [pixAcq1 spacingBetweenSlices]);
    DCMPix *pixAcq1slice2 = [pixListAcq1 objectAtIndex:1];
    printf("SliceLocation 2: %g\n", [pixAcq1slice2 sliceLocation]);
    printf("StackDirection: %i\n", [pixAcq1 stackDirection]);
    
    
    // We create a Dicom object from the source file of the slice we just selected.
    // DCMObject is the main representation of a Dicom file or attribute list for networking. The object can be an image, Structured Report, Presentation State, etc.
    // This requires the header file DCMAttributeTag.h, which I got from osirix-develop->DCM Framework. Added import line at top of DESSmapsWindowController.h
    // Setting the decodingPixelData parameter to NO simply means that the pixel data will not be decompressed or converted.
    DCMObject *dcmObjectAcq1 = [DCMObject objectWithContentsOfFile:[pixAcq1 srcFile] decodingPixelData:NO];
    
    
    // If we have only one DESS scans, we need to determine if it's a "one-touch" acquisition, i.e. two DESS scans combined into one.
    // If we have two DESS scans, then it's clearly not one-touch. So there are three possibilities:
    //
    // 1. Two separate DESS scans, with different diffusion weightings: twoSeparateDESSscans = true, twoCombinedDESSscans = false, singleDESSscan = false
    // 2. One DESS scan that is really two scans, with different diffusion weightings, combined (one-touch): twoSeparateDESSscans = false, twoCombinedDESSscans = true, singleDESSscan = false
    // 3. Just a single DESS scan, with just one diffusion weighting: twoSeparateDESSscans = false, twoCombinedDESSscans = false, singleDESSscan = true
    
    // Currently, twoCombinedDESSscan and singleDESSscan should both be false. twoSeparateDESSscans could be either true or false, depending on how many series were loaded into the plugin.
    if (twoSeparateDESSscans) {
        NSLog(@"We have two separate DESS scans.");
    } else {
        NSLog(@"We don't have two separate DESS scans.");
        // Start by assuming that it's not one-touch.
        // If the scan was one-touch, then parameters of the two DESS scans, that have been combined into one, should be stored in the Dicom header field 0043,1038. Let's read the part of the field that would contain the flip angle of the first scan, if this is one-touch.
        // This requires the header file DCMAttributeTag.h, which I got from osirix-develop->DCM Framework. Added import line at top of DESSmapsWindowController.h
        // BSv: Temporally overriding
        oneTouchAcq1flipAngle  = [[[[dcmObjectAcq1 attributeForTag:[DCMAttributeTag tagWithTagString:@"0043,1038"]] values] objectAtIndex:19] floatValue];
//        // BSv: Temporally overriding
//        oneTouchAcq1flipAngle = 25;
        if (oneTouchAcq1flipAngle > 0) {
            NSLog(@"The flip angle read from the Dicom header is %f, which is > 0, so this is one-touch.",oneTouchAcq1flipAngle);
            
            twoCombinedDESSscans = true;
        } else {
            singleDESSscan = true;
        }
    }
    

    
    // We will now read the TR, TE, and flip angle field in the Dicom header of acquisition 1.
    if ([[pixAcq1 repetitiontime] floatValue] > 0) {
        acq1tr = [[pixAcq1 repetitiontime] floatValue];
    } else {
        acq1tr = 0;
    }
    
    if ([[pixAcq1 echotime] floatValue] > 0) {
        acq1te = [[pixAcq1 echotime] floatValue];
    } else {
        acq1te = 0;
    }
    
    if ([[pixAcq1 flipAngle] floatValue] > 0) {
        acq1flipAngle = [[pixAcq1 flipAngle] floatValue];
    } else {
        acq1flipAngle = 0;
    }
    

    
    // Agx, Agy and Agz contain the areas of the gradients in each direction, units are (G/cm)*us. gradDur is the spoiler duration. For the first acquisition (or the only acquisition if there is only one), these are always stored in fields 0019,10b4-7, also for one-touch.
    float acq1Agx,acq1Agy,acq1Agz,acq1gradDur,acq1gradAmp,acq2Agx,acq2Agy,acq2Agz,acq2gradDur,acq2gradAmp;
    acq1Agx = [[[[dcmObjectAcq1 attributeForTag:[DCMAttributeTag tagWithTagString:@"0019,10b4"]] value] description] floatValue];
    acq1Agy = [[[[dcmObjectAcq1 attributeForTag:[DCMAttributeTag tagWithTagString:@"0019,10b5"]] value] description] floatValue];
    acq1Agz = [[[[dcmObjectAcq1 attributeForTag:[DCMAttributeTag tagWithTagString:@"0019,10b6"]] value] description] floatValue];
    // Duration of spoiler in ms.
    acq1gradDur = [[[[dcmObjectAcq1 attributeForTag:[DCMAttributeTag tagWithTagString:@"0019,10b7"]] value] description] floatValue]/1000;
    // Total amplitude of gradient (G/cm). First compute the total amplitude (think components of 3D triangle) and then divide by
    // duration in microseconds. Note that units of Agx/Agy/Agz are G/cm*us and units of acq1gradDur are ms. Result will be in G/cm.
    // BSv: Temporally overriding
//    acq1Agz = 1566;
//    acq1gradDur = 3.4;
    acq1gradAmp = sqrt(acq1Agx*acq1Agx + acq1Agy*acq1Agy + acq1Agz*acq1Agz)/(1000*acq1gradDur);
    

    
    if (twoSeparateDESSscans) {
        // We have another DESS scan in a separate series. We retrieve its parameters in the same way as for the 1st DESS scan.
        NSMutableArray *pixListAcq2 = [viewer pixList:2];           // Select array of slices from 2nd DESS scan
        DCMPix *pixAcq2 = [pixListAcq2 objectAtIndex:0];            // Select slice no. 0
        DCMObject *dcmObjectAcq2 = [DCMObject objectWithContentsOfFile:[pixAcq2 srcFile] decodingPixelData:NO];
        
        // Get the TR, TE, and flip angle for acquisition 2, just as we did for acquisition 1.
        if ([[pixAcq2 repetitiontime] floatValue] > 0) {
            acq2tr = [[pixAcq2 repetitiontime] floatValue];
        } else {
            acq2tr = 0;
        }
        
        if ([[pixAcq2 echotime] floatValue] > 0) {
            acq2te = [[pixAcq2 echotime] floatValue];
        } else {
            acq2te = 0;
        }
        
        if ([[pixAcq2 flipAngle] floatValue] > 0) {
            acq2flipAngle = [[pixAcq2 flipAngle] floatValue];
        } else {
            acq2flipAngle = 0;
        }
        
        // Read the gradient areas, duration, and total amplitude for acquisition 2, just as we did for acquisition 1.
        acq2Agx = [[[[dcmObjectAcq2 attributeForTag:[DCMAttributeTag tagWithTagString:@"0019,10b4"]] value] description] floatValue];
        acq2Agy = [[[[dcmObjectAcq2 attributeForTag:[DCMAttributeTag tagWithTagString:@"0019,10b5"]] value] description] floatValue];
        acq2Agz = [[[[dcmObjectAcq2 attributeForTag:[DCMAttributeTag tagWithTagString:@"0019,10b6"]] value] description] floatValue];
        // Duration in ms
        acq2gradDur = [[[[dcmObjectAcq2 attributeForTag:[DCMAttributeTag tagWithTagString:@"0019,10b7"]] value] description] floatValue]/1000;
        // Gradient amplitude in G/cm - see above.
        acq2gradAmp = sqrt(acq2Agx*acq2Agx + acq2Agy*acq2Agy + acq2Agz*acq2Agz)/(1000*acq2gradDur);
    } else if (twoCombinedDESSscans) {
        // BSv: Temporally overriding
        //acq1flipAngle = oneTouchAcq1flipAngle;
        //acq2flipAngle  = [[[[dcmObjectAcq1 attributeForTag:[DCMAttributeTag tagWithTagString:@"0043,1038"]] values] objectAtIndex:20] floatValue];
        //acq2Agx  = [[[[dcmObjectAcq1 attributeForTag:[DCMAttributeTag tagWithTagString:@"0043,1038"]] values] objectAtIndex:21] floatValue];
        //acq2Agy  = [[[[dcmObjectAcq1 attributeForTag:[DCMAttributeTag tagWithTagString:@"0043,1038"]] values] objectAtIndex:22] floatValue];
        //acq2Agz  = [[[[dcmObjectAcq1 attributeForTag:[DCMAttributeTag tagWithTagString:@"0043,1038"]] values] objectAtIndex:23] floatValue];
        
        // For one-touch scans, the TR, TE, and gradient duration should be the same for the two scans
        acq2tr = acq1tr;
        acq2te = acq1te;
        acq2gradDur = acq1gradDur;
        acq2gradAmp = sqrt(acq2Agx*acq2Agx + acq2Agy*acq2Agy + acq2Agz*acq2Agz)/(1000*acq2gradDur);
    } else {
        // This only has a single DESS acquisition (not a set of two scans, and not a "one-touch" scan). We assign values of zero to the parameters for acquisition 2.
        acq2tr = 0;
        acq2te = 0;
        acq2flipAngle = 0;
        acq2Agx = 0;
        acq2Agy = 0;
        acq2Agz = 0;
        acq2gradDur = 0;
        acq2gradAmp = 0;
    }
    
//    // BSv: Temporally overriding
//    acq1flipAngle = 25;
//    acq2flipAngle = 25;
//    acq2Agx = 0;
//    acq2Agy = 0;
//    acq2Agz = 15660;
//    acq2tr = acq1tr;
//    acq2te = acq1te;
//    acq2gradDur = acq1gradDur;
//    acq2gradAmp = sqrt(acq2Agx*acq2Agx + acq2Agy*acq2Agy + acq2Agz*acq2Agz)/(1000*acq2gradDur);
    
    
    NSLog(@"GetFieldsFromDicom has finished.");
    NSLog(@"Results:");
    NSLog(@"twoSeparateDESSscans: %@", (twoSeparateDESSscans ? @"True" : @"False"));
    NSLog(@"twoCombinedDESSscans: %@", (twoCombinedDESSscans ? @"True" : @"False"));
    NSLog(@"singleDESSscan: %@", (singleDESSscan ? @"True" : @"False"));
    NSLog(@"acq1tr: %f", acq1tr);
    NSLog(@"acq1te: %f", acq1te);
    NSLog(@"acq1flipAngle: %f", acq1flipAngle);
    NSLog(@"acq1Agx: %f", acq1Agx);
    NSLog(@"acq1Agy: %f", acq1Agy);
    NSLog(@"acq1Agz: %f", acq1Agz);
    NSLog(@"acq1gradDur: %f", acq1gradDur);
    NSLog(@"acq1gradAmp: %f", acq1gradAmp);
    NSLog(@"acq2tr: %f", acq2tr);
    NSLog(@"acq2te: %f", acq2te);
    NSLog(@"acq2flipAngle: %f", acq2flipAngle);
    NSLog(@"acq2Agx: %f", acq2Agx);
    NSLog(@"acq2Agy: %f", acq2Agy);
    NSLog(@"acq2Agz: %f", acq2Agz);
    NSLog(@"acq2gradDur: %f", acq2gradDur);
    NSLog(@"acq2gradAmp: %f", acq2gradAmp);
    
    
    
    // If we have two acquisitions, we need to find out which acquisition is the one with
    // high diffusion weighting (it should have a higher gradient amplitude).

    acq2isHigh = 0;
    if ((twoSeparateDESSscans || twoCombinedDESSscans) && (acq2gradAmp > acq1gradAmp)){
        acq2isHigh = 1;
    }
    
    //if (acq2gradAmp > acq1gradAmp){
    
    if (twoSeparateDESSscans || twoCombinedDESSscans) {         // We have two DESS scans, a "High" and "Low", presumably.
        if (acq2isHigh) {
            
            // Second acquisition has a higher diffusion weighing
            acqHiTr = acq2tr;
            acqHiTe1 = acq2te;
            acqHiTe2 = 2*acq2tr - acq2te;
            acqHiFlipAngle = acq2flipAngle;
            acqHiGradAmp = acq2gradAmp;
            acqHiGradDur = acq2gradDur;
            
            // First acquisition has a lower diffusion weighting
            acqLoTr = acq1tr;
            acqLoTe1 = acq1te;
            acqLoTe2 = 2*acq1tr - acq1te;
            acqLoFlipAngle = acq1flipAngle;
            acqLoGradAmp = acq1gradAmp;
            acqLoGradDur = acq1gradDur;
            
            // This is later used for assigning pointers to the "Hi" and "Lo" pixels in the 4D volume
            // - Here the ordering is [Lo+ Lo- Hi+ Hi-]
            //hiPixelsIndex = 2;
            //loPixelsIndex = 0;
            
        } else {
            
            // First acquisition has a higher diffusion weighing
            acqHiTr = acq1tr;
            acqHiTe1 = acq1te;
            acqHiTe2 = 2*acq1tr - acq1te;
            acqHiFlipAngle = acq1flipAngle;
            acqHiGradAmp = acq1gradAmp;
            acqHiGradDur = acq1gradDur;
            
            // Second acquisition has a lower diffusion weighting
            acqLoTr = acq2tr;
            acqLoTe1 = acq2te;
            acqLoTe2 = 2*acq2tr - acq2te;
            acqLoFlipAngle = acq2flipAngle;
            acqLoGradAmp = acq2gradAmp;
            acqLoGradDur = acq2gradDur;
            
            // This is later used for assigning pointers to the "Hi" and "Lo" pixels in the 4D volume
            // - Here the ordering is [Hi+ Hi- Lo+ Lo-]
            //hiPixelsIndex = 0;
            //loPixelsIndex = 2;
        }
    } else {        // We only have a single DESS scan. We assume it to be the "Low" scan.
        acqLoTr = acq1tr;
        acqLoTe1 = acq1te;
        acqLoTe2 = 2*acq1tr - acq1te;
        acqLoFlipAngle = acq1flipAngle;
        acqLoGradAmp = acq1gradAmp;
        acqLoGradDur = acq1gradDur;
        
        // This is later used for assigning pointers to the "Lo" pixels in the 4D volume
        // - Here the ordering is [Lo+ Lo-]
        //hiPixelsIndex = 2;      // I'm not going to use this variable, but can leave it at 2.
        //loPixelsIndex = 0;
    }
    
    NSLog(@"Have figured out which acquisition is High and which one is Low.");
    NSLog(@"acqHiTr: %f", acqHiTr);
    NSLog(@"acqHiTe1: %f", acqHiTe1);
    NSLog(@"acqHiFlipAngle: %f", acqHiFlipAngle);
    NSLog(@"acqHiGradAmp: %f", acqHiGradAmp);
    NSLog(@"acqHiGradDur: %f", acqHiGradDur);
    NSLog(@"acqLoTr: %f", acqLoTr);
    NSLog(@"acqLoTe1: %f", acqLoTe1);
    NSLog(@"acqLoFlipAngle: %f", acqLoFlipAngle);
    NSLog(@"acqLoGradAmp: %f", acqLoGradAmp);
    NSLog(@"acqLoGradDur: %f", acqLoGradDur);
}


//- (void) generateDESSmaps {
//    NSLog(@"Here I would generate a T2 map.");
//    NSLog(@"Here I would generate an ADC map.");
//}








// The closeDialog function is intended to respond when a user clicks a button.
- (IBAction) buttonPushed:(id)sender {
    //[dessWindow orderOut:sender];                               // dessWindow is a global NSWindow object, defined in DESSmapsController.h. The orderout method removes it from the screen list, hiding the window.
    //[NSApp endSheet:dessWindow returnCode:[sender tag]];          // A sheet is an instance attached to a window - see Hillegass 4th ed p. 329.
    
    //if ([sender tag] == 1) {    // User clicks OK Button
    // Do something
    
    id fittingWaitWindow = [[baseFilter viewerController] startWaitWindow:@"Producing parameter maps..."];

    float *acqHiEcho1Pix;
    float *acqHiEcho2Pix;
    float *acqLoEcho1Pix;
    float *acqLoEcho2Pix;
    float *estimatedT2slice;
    float *estimatedADCslice;
    float *testS2Pix1;
    float *testS2Pix2;
    struct diffusionArrayPoint* diffusionRatioArray;
    int sliceIndex;
    int assumedT1fromGUI_ms = (int) ([assumedT1Field intValue]);		// Units: msec
    printf("Assumed T1 from GUI (ms): %i\n", assumedT1fromGUI_ms);
    float assumedT1fromGUI = ((float) assumedT1fromGUI_ms)/1000.0;      // Units: sec
    printf("Assumed T1 from GUI (s): %g\n", assumedT1fromGUI);
    float assumedT1 = 1.2;
    int estT2type = 0;
    int estADCtype = 0;
    int slicesPerAcq;
    
    double T1start = 1.15;
    double T1end = 1.25;
    int NT1points = 2;
    double T2start = 0.01;
    double T2end = 0.08;
    int NT2points = 70;
    double Dstart = 0.1*1e-9;
    double Dend = 6*1e-9;
    int NDpoints = 60;
    
    ///////////////
//     I have temporarily given up on figuring out when and why some series get loaded with slices in the opposite order to how they are displayed in Osirix. This is important because when we are running "one-touch" scans, where the second phase slices after the first phase slices, this can cause problems. If the order gets reversed, then the plugin will actually think that the first phase is the second phase and vice versa. I did some testing for this in Exam 9270. See also Exam 33063.
//     For now, a short-term hack is to simply compute the average signal of the second echo of the two phases, in whichever order they appear. The one with the lower average S2 signal gets assigned as the high diffusion weighted ("Hi") phase, the other as the low diffusion weighted ("low") phase.
    ///////////////
    double testS2mean1 = 0;
    double testS2mean2 = 0;
    double testS2sliceSum1;
    double testS2sliceSum2;
    printf("acq2isHigh %i\n", acq2isHigh);
    int px, py;
    if (twoCombinedDESSscans){      // This is a one-touch scan - the first half of the slices corresponds to one acquisition.
        slicesPerAcq = zsize/2;
        
        printf("Start computing the sum of the two acquisitions.\n");
        for (sliceIndex = 0; sliceIndex < slicesPerAcq; sliceIndex++){
            testS2sliceSum1 = 0;
            testS2sliceSum2 = 0;
            testS2Pix1 = [[[[baseFilter viewerController] pixList:1] objectAtIndex:sliceIndex] fImage];
            testS2Pix2 = [[[[baseFilter viewerController] pixList:1] objectAtIndex:sliceIndex + slicesPerAcq] fImage];
            
            for (px=0;px<xsize;px++){
                for (py=0;py<ysize;py++){
                    testS2sliceSum1 += testS2Pix1[py*xsize+px];
                    testS2sliceSum2 += testS2Pix2[py*xsize+px];
                }
            }
            testS2mean1 = testS2sliceSum1/(xsize*ysize);
            testS2mean2 = testS2sliceSum2/(xsize*ysize);
        }
        testS2mean1 = testS2mean1/zsize;
        testS2mean2 = testS2mean2/zsize;
        printf("Done computing the mean of the two acquisitions.\n");
        printf("Mean S2 from acq1: %g\n", testS2mean1);
        printf("Mean S2 from acq2: %g\n", testS2mean2);
        
        if (testS2mean1 > testS2mean2){
            printf("acq2 seems to be High\n");
            acq2isHigh = 1;
        } else {
            printf("acq1 seems to be High\n");
            acq2isHigh = 0;
        }
    } else {
        slicesPerAcq = zsize;
    }
    
    printf("slicesPerAcq: %i\n", slicesPerAcq);
    printf("acq2isHigh %i\n", acq2isHigh);
    
    ViewerController *T2Viewer = [baseFilter duplicateCurrent2DViewerWindow];
    ViewerController *ADCViewer;
    
    
    [[[T2Viewer pixList] objectAtIndex:0] setGeneratedName:@"DESS Map T2 (ms)"];
    [[T2Viewer imageView] setWLWW:50:100];
    
    if (singleDESSscan == 0) { // We have two DESS scans, presumably a "Hi" and a "Lo" scan, so we can do an ADC fit.
        printf("Start allocating memory for the diffusion ratio array.\n");
        diffusionRatioArray = malloc(NDpoints*sizeof(struct diffusionArrayPoint));
        if (diffusionRatioArray == NULL) {
            printf("Could not allocate memory for diffusion ratio array!!!\n");
        }
        printf("Done allocating memory for the diffusion ratio array.\n");
        
        
        ADCViewer = [baseFilter duplicateCurrent2DViewerWindow];
        [[[ADCViewer pixList] objectAtIndex:0] setGeneratedName:@"DESS Map ADC (um2/s)"];
        [[ADCViewer imageView] setWLWW:1500:3000];
        
        generateDiffusionRatioArray(acqHiTr/1000, acqHiTe1/1000, acqHiFlipAngle, acqHiGradAmp*100, acqHiGradDur/1000, acqLoTr/1000, acqLoTe1/1000, acqLoFlipAngle, acqLoGradAmp*100, acqLoGradDur/1000, T1start, T1end, NT1points, T2start, T2end, NT2points, Dstart, Dend, NDpoints, diffusionRatioArray);
    }
    
    //testSlice1 = [[[[baseFilter viewerController] pixList:0] objectAtIndex:sliceIndex] fImage];

    
    
    
    for (sliceIndex = 0; sliceIndex < slicesPerAcq; sliceIndex++){
//    sliceIndex = 19;
    
        
        
        estimatedT2slice = [[[T2Viewer pixList] objectAtIndex:sliceIndex] fImage];
        
        
        if (singleDESSscan ) { // We just have one DESS scan - assume that it's the "Lo" scan, with low spoiling.
            acqLoEcho1Pix = [[[[baseFilter viewerController] pixList:0] objectAtIndex:sliceIndex] fImage];
            acqLoEcho2Pix = [[[[baseFilter viewerController] pixList:1] objectAtIndex:sliceIndex] fImage];
            NSLog(@"Source image pointers have been set for single DESS scan.");
        } else {
            if (twoCombinedDESSscans == 0){  // This is not a "one-touch" scan
                if (acq2isHigh) {
                    acqLoEcho1Pix = [[[[baseFilter viewerController] pixList:0] objectAtIndex:sliceIndex] fImage];
                    acqLoEcho2Pix = [[[[baseFilter viewerController] pixList:1] objectAtIndex:sliceIndex] fImage];
                    acqHiEcho1Pix = [[[[baseFilter viewerController] pixList:2] objectAtIndex:sliceIndex] fImage];
                    acqHiEcho2Pix = [[[[baseFilter viewerController] pixList:3] objectAtIndex:sliceIndex] fImage];
                    NSLog(@"Source image pointers have been set for separate Lo-Hi DESS scans (in that order).");
                } else {
                    acqHiEcho1Pix = [[[[baseFilter viewerController] pixList:0] objectAtIndex:sliceIndex] fImage];
                    acqHiEcho2Pix = [[[[baseFilter viewerController] pixList:1] objectAtIndex:sliceIndex] fImage];
                    acqLoEcho1Pix = [[[[baseFilter viewerController] pixList:2] objectAtIndex:sliceIndex] fImage];
                    acqLoEcho2Pix = [[[[baseFilter viewerController] pixList:3] objectAtIndex:sliceIndex] fImage];
                    NSLog(@"Source image pointers have been set for separate Hi-Lo DESS scans (in that order).");
                }
            } else { // This is a one-touch scan
                
                
                if (acq2isHigh) {
//                if (acq2isHigh == 0) {
                    acqLoEcho1Pix = [[[[baseFilter viewerController] pixList:0] objectAtIndex:sliceIndex] fImage];
                    acqLoEcho2Pix = [[[[baseFilter viewerController] pixList:1] objectAtIndex:sliceIndex] fImage];
                    acqHiEcho1Pix = [[[[baseFilter viewerController] pixList:0] objectAtIndex:sliceIndex + slicesPerAcq] fImage];
                    acqHiEcho2Pix = [[[[baseFilter viewerController] pixList:1] objectAtIndex:sliceIndex + slicesPerAcq] fImage];
                    NSLog(@"Source image pointers have been set for one-touch Lo-Hi DESS scans (in that order).");
                } else {
                    acqHiEcho1Pix = [[[[baseFilter viewerController] pixList:0] objectAtIndex:sliceIndex] fImage];
                    acqHiEcho2Pix = [[[[baseFilter viewerController] pixList:1] objectAtIndex:sliceIndex] fImage];
                    acqLoEcho1Pix = [[[[baseFilter viewerController] pixList:0] objectAtIndex:sliceIndex + slicesPerAcq] fImage];
                    acqLoEcho2Pix = [[[[baseFilter viewerController] pixList:1] objectAtIndex:sliceIndex + slicesPerAcq] fImage];
                    NSLog(@"Source image pointers have been set for one-touch Hi-Lo DESS scans (in that order).");
                }
            }
            
            //(364,243)
//            printf("sliceIndex: %i\n", sliceIndex);
//            printf("acqHiEcho1Pix[377,74]: %g\n", acqHiEcho1Pix[77*512+374]);
//            printf("acqHiEcho2Pix[377,74]: %g\n", acqHiEcho2Pix[77*512+374]);
//            printf("acqLoEcho1Pix[377,74]: %g\n", acqLoEcho1Pix[77*512+374]);
//            printf("acqLoEcho2Pix[377,74]: %g\n", acqLoEcho2Pix[77*512+374]);
//
            
//            int xx, yy;
//            for (xx=100;xx<200;xx++){
//                for (yy=100;yy<200;yy++){
//                    acqHiEcho1Pix[yy*xsize+xx] = 1234;
//                    acqHiEcho2Pix[yy*xsize+xx] = 1234;
//                    acqLoEcho1Pix[(yy+50)*xsize+(xx+50)] = 1234;
//                    acqLoEcho2Pix[(yy+50)*xsize+(xx+50)] = 1234;
//                }
//            }
            
        }
    
//        // Find the maximum value of the S1L signal, to set a noise threshold (something like 0.15*S1L)
//        int pp;
//        float maxS1L = 0;
//        for (pp = 0; pp < xsize*ysize; pp++){
//            if (acqLoEcho1Pix[pp] > maxS1L){
//                maxS1L = acqLoEcho1Pix[pp];
//            }
//        }
        
        
//        computeT2map(acqLoTr/1000, acqLoTe1/1000, acqLoFlipAngle, assumedT1, acqLoEcho1Pix, acqLoEcho2Pix, xsize*ysize, estT2type, 0.15*maxS1L, estimatedT2slice);
//        computeT2map(acqLoTr/1000, acqLoTe1/1000, acqLoFlipAngle, assumedT1, acqLoEcho1Pix, acqLoEcho2Pix, xsize*ysize, estT2type, estimatedT2slice);
        computeT2map(acqLoTr/1000, acqLoTe1/1000, acqLoFlipAngle, assumedT1fromGUI, acqLoEcho1Pix, acqLoEcho2Pix, xsize*ysize, estT2type, estimatedT2slice);
        NSLog(@"Done computing T2 map.");
        maskNoise(estimatedT2slice, acqLoEcho1Pix, 0.15, xsize, ysize);
//        printf("estimatedT2[95,135]: %g\n", estimatedT2slice[135*256+95]*1000);
        //computeT2map(float TR, float TE, float alphaL, float assumedT1, float* SL1, float* SL2, int imageSize, int estType, float*    estimatedT2)
        
        if (singleDESSscan == 0){ // We have two DESS scans, presumably a "Hi" and a "Lo" scan, so we can do an ADC fit.
            estimatedADCslice = [[[ADCViewer pixList] objectAtIndex:sliceIndex] fImage];
//            double ReS1H,ImS1H,ReS2H,ImS2H,ReS1L,ImS1L,ReS2L,ImS2L;
//            ReS1H = 0;
//            ImS1H = 0;
//            ReS2H = 0;
//            ImS2H = 0;
//            ReS1L = 0;
//            ImS1L = 0;
//            ReS2L = 0;
//            ImS2L = 0;
//            computeEchoesEPG(1.2,0.04,0.018,0.005,30,400,0.0034,1.5*1e-9,6,&ReS1H,&ImS1H,&ReS2H,&ImS2H);
//            printf("S1H: %g+i%.15g\n", ReS1H, ImS1H);
//            printf("S2H: %g+i%.15g\n", ReS2H, ImS2H);
//            computeEchoesEPG(1.2,0.04,0.018,0.005,30,4,0.0034,1.5*1e-9,6,&ReS1L,&ImS1L,&ReS2L,&ImS2L);
//            printf("S1H: %g+i%.15g\n", ReS1L, ImS1L);
//            printf("S2H: %g+i%.15g\n", ReS2L, ImS2L);
//            printf("SCH: %g\n", SC(ReS2H,ImS2H,ReS1H,ImS1H));
//            printf("SCL: %g\n", SC(ReS1L,ImS1L,ReS2L,ImS2L));
//            printf("SCH*SCL: %g\n", SC(ReS2H,ImS2H,ReS1H,ImS1H)*SC(ReS1L,ImS1L,ReS2L,ImS2L));
//            computeADCmap(acqHiTr/1000, acqHiTe1/1000, acqHiFlipAngle, acqHiGradAmp*100, acqHiGradDur/1000, acqLoTr/1000, acqLoTe1/1000, acqLoFlipAngle, acqLoGradAmp*100, acqLoGradDur/1000, acqHiEcho1Pix, acqHiEcho2Pix, acqLoEcho1Pix, acqLoEcho2Pix, xsize*ysize, estADCtype, 0.15*maxS1L, estimatedADCslice);
            computeADCmap(acqHiEcho1Pix, acqHiEcho2Pix, acqLoEcho1Pix, acqLoEcho2Pix, xsize*ysize, estADCtype,estimatedADCslice, diffusionRatioArray, NDpoints);
            NSLog(@"Done computing ADC map.");
            maskNoise(estimatedADCslice, acqLoEcho1Pix, 0.15, xsize, ysize);

        }
    }
    [[baseFilter viewerController] endWaitWindow: fittingWaitWindow];
    NSLog(@"Estimate button was pushed");
    //[self generateDESSmaps];
    //}
}





// The dealloc function is standard for deallocating memory. Must release any objects that the object was retaining and the call the superclass's dealloc method (Hillegass, 4th ed, p. 72).
- (void) dealloc {
    [super dealloc];
}





// The windowWillClose function prevents memory leaks, by releasing memory after the window is closed (I think).
- (void) windowWillClose:(NSNotification*)notification {
    if ([notification object] == [self window]) {
        [[self window] orderOut: self];
        [self autorelease];
    }
}



@end
