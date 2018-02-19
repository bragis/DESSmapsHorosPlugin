//
//  DESSmapsFilter.h
//  DESSmaps
//
//  Copyright (c) 2016 Bragi. All rights reserved.
//


#import <Foundation/Foundation.h>
#import <OsiriXAPI/PluginFilter.h>

@interface DESSmapsFilter : PluginFilter {

}

- (long) filterImage:(NSString*) menuName;
//- (IBAction) closeDialog:(id)sender;
- (ViewerController*) viewerController;

@end