//
//  DESSmapsFilter.m
//  DESSmaps
//
//  Copyright (c) 2016 Bragi. All rights reserved.
//

#import "DESSmapsFilter.h"
#import "DESSmapsWindowController.h"

@implementation DESSmapsFilter

- (void) initPlugin
{
}

- (ViewerController*)   viewerController
{
    return viewerController;
}

//// The closeDialog function is intended to respond when a user clicks a button.
//- (IBAction) closeDialog:(id)sender {
//    [window orderOut:sender];                               // window is a global NSWindow object, defined in DESSmapsFilter.h
//    [NSApp endSheet:window returnCode:[sender tag]];
//    
//    if ([sender tag] == 1) {    // User clicks OK Button
//        // Do something
//    }
//}

- (long) filterImage:(NSString*) menuName
{
    // We could initialize the DESSmapsWindowController like this:
    //DESSmapsWindowController *c = [[DESSmapsWindowController alloc] initWithWindowNibName: @"DESSmapsWindow"];
    //However, we want to use our own initalization method, initWithFilter, that calls initWithWindowNibName and does some other stuff too. So we use this line:
    
    DESSmapsWindowController	*c = [[DESSmapsWindowController alloc] initWithFilter:self];
    
    //I'm not sure what the difference is between the line below and "[c showWindow:self]", but the below seems to work.
    [[c window] makeKeyAndOrderFront: self];
    
    return 0;
}

@end
