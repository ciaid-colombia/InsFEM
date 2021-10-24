#!/usr/bin/perl -w
use FileHandle;
use Data::Dumper;

my $postfile = "Tenerife_caso34_v1.post.res";
my $postfolder = "figures";

#Time stepping
my $TotalNumberOfSteps = 500;
my $PostprocessEvery   = 10;
#Video
my %Video;
$Video{Do} = "Yes";
$Video{Format} = "AVIMSVC";
$Video{FramesPerStep} = 5;

#Fields
my @Fields;

$Fields[0]{Name} = "TEMPE";
$Fields[0]{IsLog}="No";
$Fields[0]{Limits}{IsLimited}="No";$Fields[0]{Limits}{0} = 24.5; $Fields[0]{Limits}{1} = 34.1;

$Fields[1]{Name} = "VELOC |VELOC|";

$Fields[2]{Name} = "CN2";
$Fields[2]{IsLog}="Yes";

$Fields[3]{Name} = "PRESS";

my @Volumes;
$Volumes[0] = "V_Tetrahedra";

#Cuts
my @Cuts;

$Cuts[0]{Name} = "C_CutSet_1_Tetrahedra";
$Cuts[0]{Coords} = "0 0 2428 1 0 2428 0 1 2428";
$Cuts[0]{View} = "Topview.vv";


$Cuts[1]{Name} = "C_CutSet_2_Tetrahedra";
$Cuts[1]{Coords} = "0 0 0 0 1 0 0 0 1";
$Cuts[1]{View} = "cut3.vv";

$Cuts[2]{Name} = "C_CutSet_3_Tetrahedra";
$Cuts[2]{Coords} = "0 0 0 1 0 0 0 0 1";
$Cuts[2]{View} = "cut2.vv";







#We start postprocessing
system ("mkdir -p $postfolder");

#First we write the file

open(OUT,">","./GidBchFile.bch");
print OUT <<EOBCH;
Mescape Postprocess
Mescape Files Read $postfile
Mescape DisplayStyle bodybound
EOBCH
#Do Cuts
for ($iCut = 0; $iCut < @Cuts; $iCut = $iCut +1) {
   print OUT "Mescape DoCut CutPlane ThreePoints $Cuts[$iCut]{Coords} escape\n";
}
#UnSelect All cuts and volumes
for my $volu (@Volumes) {
   print OUT "Select VolumeSets $volu escape escape escape escape escape\n";
   
}
for ($iCut = 0; $iCut < @Cuts; $iCut = $iCut +1) {
   print OUT "Select CutSets $Cuts[$iCut]{Name} escape escape escape escape escape escape\n";
}

#Go to the last step
#For each field, take a snapshot at the last step, forall views and cuts
for ($iField = 0; $iField < @Fields; $iField = $iField+1) {
   print OUT "Mescape Results ContourFill $Fields[$iField]{Name}\n";
   if ($Fields[$iField]{IsLog} eq "Yes") {
      print OUT "Mescape Results Options ScaleResult UseLogScale Yes escape\n";
   }
   if ($Fields[$iField]{Limits}{IsLimited} eq "Yes") {
      print OUT "Mescape results contoptions setminoptions setvalue $Fields[$iField]{Limits}{0}  Mescape\n";
      print OUT "Mescape results contoptions setmaxoptions setvalue $Fields[$iField]{Limits}{1} Mescape\n";
      print OUT "Mescape Utilities Redisplay escape escape\n";
   }
   
  
   print OUT "Mescape Results ContourFill $Fields[$iField]{Name}\n";
   #Loop through cuts
   for ($iCut = 0; $iCut < @Cuts; $iCut = $iCut +1) {
      print OUT "Mescape results analysissel ANALYSIS $TotalNumberOfSteps escape\n";
      print OUT "escape escape escape escape\n";
      print OUT "Select CutSets $Cuts[$iCut]{Name} escape escape escape escape escape\n";
      print OUT "escape escape escape escape\n";
      print OUT "Mescape escape escape escape escape escape escape View ReadView $Cuts[$iCut]{View}\n";
      print OUT "\'HardCopy PNG \"$postfolder/$Fields[$iField]{Name}_Cut$iCut.png\" \n";
      print OUT "escape escape escape escape\n";
      
      #Videos
      if ($Video{Do} eq "Yes") {
         print OUT <<EOBCH;
escape escape escape escape Utilities Variables PostUpdateWindows yes escape escape escape escape
'AnimationFile Format $Video{Format}
'AnimationFile FramesPerStep $Video{FramesPerStep} escape escape escape escape
'AnimationFile Start "$postfolder/$Fields[$iField]{Name}_Cut$iCut.avi" Yes escape escape escape escape
escape escape escape escape Utilities Variables PostUpdateWindows no escape escape escape escape
EOBCH
         for ($istep = $PostprocessEvery; $istep <= $TotalNumberOfSteps; $istep = $istep + $PostprocessEvery) {
            print OUT "Mescape Results AnalysisSel ANALYSIS $istep Mescape\n";
            print OUT "Mescape Results ContourFill $Fields[$iField]{Name}\n";              
            print OUT "\'AnimationFile AddStep\n";
         }
         print OUT "escape escape escape escape Utilities Variables PostUpdateWindows yes escape escape escape escape\n";
         print OUT "\'AnimationFile End\n";

      }
      
      #Unselect cut
      print OUT "Select CutSets $Cuts[$iCut]{Name} escape escape escape escape escape\n";
      print OUT "escape escape escape escape\n";
      

   }
   
   
   print OUT "Mescape Results Options ScaleResult UseLogScale NO escape\n";
}
print OUT "quit\n";

close OUT;

system("gid -b+g GidBchFile.bch");
unlink "GidBchFile.bch";






