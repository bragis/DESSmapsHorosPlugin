//
//  DESSmapsWindowController.h
//  DESSmaps
//
//  Created by Bragi Sveinsson on 1/16/16.
//
//

#import <Cocoa/Cocoa.h>
#import "DESSmapsFilter.h"
#import "DCMObject.h"
#import "DCMAttributeTag.h"
#import "estimateT2.h"
#import "estimateADC.h"
#import "noiseMasking.h"

@interface DESSmapsWindowController : NSWindowController {
    DESSmapsFilter *baseFilter;         // Viewer that called the plugin.
    
    // Image Sizes (get/calculate for convenience)
    long nxyzpts;	// total # of DESS pixels
    long xsize;		// # of DESS pixels in X
    long ysize;		// # of DESS pixels in Y
    long zsize;		// # of DESS pixels (slices) in Z
    long ntpts;		// # of DESS time points (or series in 4D viewer).
    
    bool twoSeparateDESSscans;      // Used to be twoDESSscans
    bool twoCombinedDESSscans;      // Used to be isOneTouch
    bool singleDESSscan;
    bool acq2isHigh;
    
    float acqHiTr;
    float acqHiTe1;
    float acqHiTe2;
    float acqHiFlipAngle;
    float acqHiGradAmp;
    float acqHiGradDur;
    float acqLoTr;
    float acqLoTe1;
    float acqLoTe2;
    float acqLoFlipAngle;
    float acqLoGradAmp;
    float acqLoGradDur;
    
    //int hiPixelsIndex;
    //int loPixelsIndex;
    
    /*acqHiTr = acq1tr;
    acqHiTe1 = acq1te;
    acqHiTe2 = 2*acq1tr - acq1te;
    acqHiFlipAngle = acq1flipAngle;
    acqHiGradAmp = acq1gradAmp;
    acqHiGradDur = acq1gradDur; */
    
    IBOutlet NSTextField *acqHiTrField;
    IBOutlet NSTextField *acqHiTe1Field;
    IBOutlet NSTextField *acqHiTe2Field;
    IBOutlet NSTextField *acqFlipAngleField;
    IBOutlet NSTextField *acqGradAmpField;
    IBOutlet NSTextField *acqGradDurField;
    IBOutlet NSTextField *assumedT1Field;
    
    IBOutlet NSWindow *dessWindow;      // The main window of the app.
}

- (id) initWithFilter:(DESSmapsFilter *) filter;
- (void) getFieldsFromDicom;

- (void)updateFieldsInGUI;

- (IBAction) buttonPushed:(id)sender;

@end
