use strict;
use MdmDiscoveryScript;
use SdmDiscoveryScript;
use ProtocolCommands;
use ProtocolDiscoveryScript;
use File::Path qw(rmtree);
use ForceFieldDiscoveryScript;
use File::Basename;
use warnings;
use Cwd;
use File::Copy 'move';
use ChartDiscoveryScript;
use DSScript;
use ProteinDiscoveryScript;
use utf8;
print cwd()."\n";
my $work_dir = 'E:\00-surf-anaylsis-HPB-Chrg\struct';
$work_dir =~ s/\\/\//g;
printf "$work_dir\n";

sub sub_run_calculate_CDR_IMGT_HPB{ 
    my $filename = $_[0];
    my $work_dir = $_[1];
    my @suffixlist = qw(.pdb .pl .txt);
    my $basename = basename($filename, @suffixlist);
    print "$basename\n";

    # Open MDM document in Molecule Window
    my $mdmDocument = DiscoveryScript::Open(
    {
    Path => $filename,
    ModelType => MdmModelType
    }
    );

    printf "Document contains Molecules: %d Total Atoms: %d\n", $mdmDocument->Molecules->Count,$mdmDocument->Atoms->Count;
    sleep(1);

    #By default, protein sequences are extracted
    my $sdmDocument = Sdm::Document::CreateDocumentFromStructures($mdmDocument);
    my $count = $sdmDocument->{SequenceCount};
    print "First document contains $count protein sequences.\n";
    my $sequencesAlignment = $sdmDocument->Alignment;


#01 ######################### START CDR annotation protocol ,setup parameter and running by script #######################
    # Protocol name and parameters
    my $protocolName = 'Annotate Antibody Sequence';
    my %parameterHash = ( "Input Sequences" =>  $sequencesAlignment,
                        "Annotation Scheme" => "IMGT",  
                        "Renumber Sequences" =>False, 
                        );

    # print "Connecting to server...\n";
    my $server = "localhost:9943";
    print "Connecting to server $server...\n";
    my $session = Protocol::Session::Create($server);
    # my $session = Protocol::Document::DefaultSession();
    if ($session)
    {
        print "Connected to $session->{Server} as user $session->{User}.\n";
    }
    else
    {
        die "No Pipeline Pilot server found. Please connect to a valid server and try again.\n";
    }

    # Print default running folder, which can be changed
    print "Running in $session->{RunFolder}...\n";

    print "Open protocol document for $protocolName\n";
    my $document = Protocol::Document::Create( $protocolName, $session );
    
    my $guid = $document->Guid();
    print "Protocol guid: $guid\n"; 

    my $version = $document->Version();
    print "Protocol version: $version\n";
    ## Setup the protocol parameters

    for my $key (keys %parameterHash) {
        if ( $document->ItemExists($key) ) {
            $document->ReplaceItem( $key, $parameterHash{$key} );
            print "Set parameter '$key' : $parameterHash{$key}\n"; 
        }
        else{
            print "\n<warning> not a valid parameter '$key'\n\n"; 
        }
    }

    ## Launch Pipeline Pilot protocol

    my $saveIntermediateFiles = True;
    my $maxFilesizeToDownload = 256;
    my $showJobCompleteDialog = False;
    my $task                  = $session->Launch(
        $document,              $maxFilesizeToDownload,
        $saveIntermediateFiles, $showJobCompleteDialog
    );

    print "\nLaunch the protocol: $protocolName \n";
    print "Task info\n";
    print "Name:    $task->{ProtocolName}\n";
    print "Host:    $task->{Host}\n";
    print "RunPath: $task->{RunPath}\n";
    print "TaskId:  $task->{TaskId}\n";

    ## Track progress of task
    print " $protocolName Task is running...\n";
    my $times = 0;
    while ( $task->IsRunnable() )
    {
        print "Status: $task->{State}\n";
        sleep(5);
        ++$times;
        if ( $times > 5 )
        {
            $task->KillTask();
        }
    }

    $task->WaitForCompletion();
#######################################  END of CDR annotation ###########################################

#02 ###################################### Start cdr selection ###############################################
    my $upScheme = "IMGT";
    my $scheme = lc($upScheme);
    my $msg = 
    "Given antibody sequences annotated using the $upScheme annotation scheme, this \nscript will update the coloring in the sequence window to match this annotation.\n";

    my $path_cdrannotation = $task->RunPath;
    my $outFile = "$path_cdrannotation/Output/Sequence.bsml";
    die "Output file was not found please check report!" if ( not -e $outFile );
    # load the results        
    # DiscoveryScript::Open( { Path => $outFile, LoadAllObjects => True } ); 
    my $cdrDocument = DiscoveryScript::Open( { Path => $outFile, LoadAllObjects => True } ); 

    die "$msg\nNo open sequence window was found. \n" unless $cdrDocument;
    my $sequences = $cdrDocument->Sequences;
    print "selected seqs Counts:".$sequences->Count ."\n";
    # print $sequences->Name;
    my $schemeFound = False;
    # my $aminoAcids = Mdm::Array::Create();

    my @cdr_chain_selection = ();
    my @cdr_residue_selection = ();

    foreach my $sequence (@$sequences) {
        my $seqname = uc($sequence->Name);
        print "seqname:".$seqname ."\n";
        unless ($seqname =~ "ALIGNMENT_CONSENSUS") {
            my $count = $sequence->FeatureCount;
            print "feature counts:". $count."\n";

            my $cdrNewArray = Sdm::Array::Create();
            
            for ( my $index = 0 ; $index < $count ; ++$index ) {
                my $feature  = $sequence->Feature($index);
                my $key = $feature->Key;
                print "feature key:". $key."\n";
    #print "feature key " . $key;
                if ($key =~ "cdr" and $key =~ $scheme) {
                    $cdrNewArray->AddItem($feature);
                }
            }


            if($cdrNewArray->Count >0){
                $schemeFound = True;
                print "Changing colors on the following features for ". $sequence->Name .":\n";
        
    # add cdr residues to array of aminoacds for the follwing loop 
            foreach my $usrFeat (@$cdrNewArray) {
                print $usrFeat->Key."\n";
                my $location = $usrFeat->Location;  
                my $segment = $location->Segment(0);
                my $start = $segment->StartPosition;
                my $end = $segment->EndPosition;

                push(@cdr_chain_selection,$usrFeat->Key);
                push(@cdr_residue_selection,$sequence->Residue($start)->{Id}."-".$sequence->Residue($end)->{Id});
                }
        }
        }   

    }
    # $document->UpdateViews();
    print "@cdr_chain_selection\n";
    print "@cdr_residue_selection\n";
#################################### END test cdr selection ###############################################

#03 ################################### START of HPB calculation of IMGT CDR ##############################################

    # added 20260418 :select IMGT CDRS;Create an array of all aminoacids in the system
    $mdmDocument->SelectByProperty('Chain name = "H";Residue id="'.$cdr_residue_selection[0].'",id="'.$cdr_residue_selection[1].'",id="'.$cdr_residue_selection[2].'"');#正确语法案例
    $mdmDocument->SelectByProperty('Chain name = "L";Residue id="'.$cdr_residue_selection[3].'",id="'.$cdr_residue_selection[4].'",id="'.$cdr_residue_selection[5].'"');#正确语法案例
    # $mdmDocument->SelectByProperty('Chain name= "H";Residue id="25-32",id="50-56",id="95-104"');#正确语法案例
    my $aminoAcids = $mdmDocument->Filter({Selected => 1 })->AminoAcids;

    # Calculate the solvent accessible surface per residue
    $mdmDocument->SolventAccessibleSurface({GridNumber => 240, ContextSelected => False});

    # Loop over array of aminoacids.
    my $totalSAS = 0.0;
    my $totalHPBSAS = 0.0;
    my $HPB_residues = Mdm::Array::Create();
    my $totalY = 0.0;
    my $totalW = 0.0;
    my $total_hpb_Y_W = 0.0;
    my $YW_residues = Mdm::Array::Create();
    my $HPB_AND_YW_residues = Mdm::Array::Create();
    foreach my $aminoacid (@$aminoAcids)
    {
        # retrieve the residue name
        my $aminoAcidName = $aminoacid->Name;

        # retrieve the Residue Solvent Accessibility
        my $scSAS = $aminoacid->GetProperty('Sidechain Solvent Accessibility');
        my $SAS = $aminoacid->GetProperty('Residue Solvent Accessibility');
        $totalSAS += $SAS;
        #    print "$aminoAcidName : $scSAS $SAS\n";
        
        my $hpb = $aminoacid->GetProperty('Hydrophobicity');
            
        if ( $hpb > 0.0 )
        {
            $totalHPBSAS += $scSAS;
            if ($scSAS > 5.0) { $HPB_residues->AddItem($aminoacid);$HPB_AND_YW_residues->AddItem($aminoacid) }
        }
        
        if ( grep /^TYR|^TRP/, $aminoAcidName )
        {
            if ( $aminoacid->GetProperty('PDB Name') eq 'TYR')
            { $totalY += $scSAS;}
            if ( $aminoacid->GetProperty('PDB Name') eq 'TRP')
            { $totalW += $scSAS; }
            if ($scSAS > 5.0) { $YW_residues->AddItem($aminoacid); $HPB_AND_YW_residues->AddItem($aminoacid)}
        }        
    }

    $total_hpb_Y_W = $totalHPBSAS + $totalY + $totalW;
    $total_hpb_Y_W = sprintf("%.2f",$total_hpb_Y_W);
    $mdmDocument->CreateGroup("CDR IMGT HPB+Y+W Sidechain: $total_hpb_Y_W;", $HPB_AND_YW_residues);

    $totalSAS = sprintf("%.2f", $totalSAS);
    $totalHPBSAS = sprintf("%.2f", $totalHPBSAS);
    print "Total Solvent Accessibility : $totalSAS\n";
    print "Total Hydrophobic Sidechain Solvent Accessibility : $totalHPBSAS\n";
    $mdmDocument->CreateGroup("CDR IMGT HPB Sidechain: $totalHPBSAS; Total: $totalSAS", $HPB_residues);

    $totalY = sprintf("%.2f", $totalY);
    $totalW = sprintf("%.2f", $totalW);
    print "Total TYR Solvent Accessibility : $totalY\n";
    print "Total TRP Solvent Accessibility : $totalW\n";
    $mdmDocument->CreateGroup("CDR IMGT Y Sidechain: $totalY; CDR IMGT W Sidechain: $totalW", $YW_residues);



    # remove file
    # unlink ($path."/*");
    # rmtree($path_cdrannotation, { verbose => 0, keep_root => 0 });
#################################### END of HPB calculation of IMGT CDR ##############################################

#04 ################################### START of HPB calculation of full molecule ##############################################
    $mdmDocument->SelectAll();
    my $v2_aminoAcids = $mdmDocument->Filter({Selected => 1 })->AminoAcids;

    # Calculate the solvent accessible surface per residue
    $mdmDocument->SolventAccessibleSurface({GridNumber => 240, ContextSelected => False});

    # Loop over array of aminoacids.
    my $v2_totalSAS = 0.0;
    my $v2_totalHPBSAS = 0.0;
    my $v2_HPB_residues = Mdm::Array::Create();
    my $v2_totalY = 0.0;
    my $v2_totalW = 0.0;
    my $v2_total_hpb_Y_W = 0.0;
    my $v2_YW_residues = Mdm::Array::Create();
    my $v2_HPB_AND_YW_residues = Mdm::Array::Create();
    foreach my $v2_aminoacid (@$v2_aminoAcids)
    {
        # retrieve the residue name
        my $v2_aminoAcidName = $v2_aminoacid->Name;

        # retrieve the Residue Solvent Accessibility
        my $v2_scSAS = $v2_aminoacid->GetProperty('Sidechain Solvent Accessibility');
        my $v2_SAS = $v2_aminoacid->GetProperty('Residue Solvent Accessibility');
        $v2_totalSAS += $v2_SAS;
        #    print "$v2_aminoAcidName : $v2_scSAS $v2_SAS\n";
        
        my $v2_hpb = $v2_aminoacid->GetProperty('Hydrophobicity');
            
        if ( $v2_hpb > 0.0 )
        {
            $v2_totalHPBSAS += $v2_scSAS;
            if ($v2_scSAS > 5.0) { $v2_HPB_residues->AddItem($v2_aminoacid);$v2_HPB_AND_YW_residues->AddItem($v2_aminoacid) }
        }
        
        if ( grep /^TYR|^TRP/, $v2_aminoAcidName )
        {
            if ( $v2_aminoacid->GetProperty('PDB Name') eq 'TYR')
            { $v2_totalY += $v2_scSAS;}
            if ( $v2_aminoacid->GetProperty('PDB Name') eq 'TRP')
            { $v2_totalW += $v2_scSAS; }
            if ($v2_scSAS > 5.0) { $v2_YW_residues->AddItem($v2_aminoacid); $v2_HPB_AND_YW_residues->AddItem($v2_aminoacid)}
        }        
    }

    $v2_total_hpb_Y_W = $v2_totalHPBSAS + $v2_totalY + $v2_totalW;
    $v2_total_hpb_Y_W = sprintf("%.2f",$v2_total_hpb_Y_W);
    $mdmDocument->CreateGroup("full Mol HPB+Y+W Sidechain: $v2_total_hpb_Y_W;", $v2_HPB_AND_YW_residues);

    $v2_totalSAS = sprintf("%.2f", $v2_totalSAS);
    $v2_totalHPBSAS = sprintf("%.2f", $v2_totalHPBSAS);
    print "Total Solvent Accessibility : $v2_totalSAS\n";
    print "Total Hydrophobic Sidechain Solvent Accessibility : $v2_totalHPBSAS\n";
    $mdmDocument->CreateGroup("full Mol Hydrophobic Sidechain: $v2_totalHPBSAS; Total: $v2_totalSAS", $v2_HPB_residues);

    $v2_totalY = sprintf("%.2f", $v2_totalY);
    $v2_totalW = sprintf("%.2f", $v2_totalW);
    print "Total TYR Solvent Accessibility : $v2_totalY\n";
    print "Total TRP Solvent Accessibility : $v2_totalW\n";
    $mdmDocument->CreateGroup("full Mol Y Sidechain: $v2_totalY; W Sidechain: $v2_totalW", $v2_YW_residues);

#################################### END of HPB calculation of full molecule  ##############################################

#05 ################################### START aggregation score calculate  ###############################################
    # Protocol name and parameters
    my $v3_molecule = $mdmDocument ->Molecules->Item(0);

    print "Creating CHARMm force field document.\n";
    my $forcefield = Ffdm::Document::Create(Ffdm::charmm36ForceField);

    print "Applying CHARMm force field to mdmDocument\n";
    my $status = $forcefield->ApplyForceField( $mdmDocument, True, True, True,
        Ffdm::useMomanyRonePartialCharges );

    if ( $status == Ffdm::applyForceFieldCompleted )
    {
        print "Molecule successfully typed.\n";
    }
    else
    {
        die "Failed to type molecule.\n";
    }

    my $v3_protocolName = 'Calculate Aggregation Scores';
    my %v3_parameterHash = ( "Input Typed Protein" => $v3_molecule,
                        "Cutoff Radius" => "5,10",  
                        );

    # print "Connecting to server...\n";
    my $v3_server = "localhost:9943";
    print "Connecting to server $v3_server...\n";
    my $v3_session = Protocol::Session::Create($v3_server);
    # my $v3_session = Protocol::Document::DefaultSession();
    if ($v3_session)
    {
        print "Connected to $v3_session->{Server} as user $v3_session->{User}.\n";
    }
    else
    {
        die "No Pipeline Pilot server found. Please connect to a valid server and try again.\n";
    }

    # Print default running folder, which can be changed
    print "Running in $v3_session->{RunFolder}...\n";

    print "Open protocol document for $v3_protocolName\n";
    my $v3_document = Protocol::Document::Create( $v3_protocolName, $v3_session );
    
    my $v3_guid = $v3_document->Guid();
    print "Protocol guid: $v3_guid\n"; 

    my $v3_version = $v3_document->Version();
    print "Protocol version: $v3_version\n";
    ## Setup the protocol parameters

    for my $v3_key (keys %v3_parameterHash) {
        if ( $v3_document->ItemExists($v3_key) ) {
            $v3_document->ReplaceItem( $v3_key, $v3_parameterHash{$v3_key} );
            print "Set parameter '$v3_key' : $v3_parameterHash{$v3_key}\n"; 
        }
        else{
            print "\n<warning> not a valid parameter '$v3_key'\n\n"; 
        }
    }

    ## Launch Pipeline Pilot protocol

    my $v3_saveIntermediateFiles = True;
    my $v3_maxFilesizeToDownload = 256;
    my $v3_showJobCompleteDialog = False;
    my $v3_task                  = $v3_session->Launch(
        $v3_document,              $v3_maxFilesizeToDownload,
        $v3_saveIntermediateFiles, $v3_showJobCompleteDialog
    );

    print "\nLaunch the protocol: $v3_protocolName \n";
    print "Task info\n";
    print "Name:    $v3_task->{ProtocolName}\n";
    print "Host:    $v3_task->{Host}\n";
    print "RunPath: $v3_task->{RunPath}\n";
    print "TaskId:  $v3_task->{TaskId}\n";

    ## Track progress of task
    print " $v3_protocolName Task is running...\n";
    my $v3_times = 0;
    while ( $task->IsRunnable() )
    {
        print "Status: $task->{State}\n";
        sleep(5);
        ++$times;
        if ( $times > 5 )
        {
            $task->KillTask();
        }
    }

    $v3_task->WaitForCompletion();
    $mdmDocument->UpdateViews();
    $sdmDocument->UpdateViews();
    # remove file
    # unlink ($path."/*");
   
#################################### END aggregation score calculate  ###############################################

#06 ################################ START aggregation score visulization  ###############################################
    my $path_agg_calculation = $v3_task->RunPath;
    $mdmDocument->Insert({"Path" => "$path_agg_calculation/Output/$basename.mol2", "FormatType" => "mol2-merge-properties"});
    #make all but this molecule invisible and make it the defined reference
    #this only works if OBJECTID is defined.
    my $molecules = $mdmDocument->Molecules;
    my $nmol = $mdmDocument->Molecules->Count;
    print "Mol count:$nmol";
    my $refMol = -1;
    #want to make sure we get strings to compare.
    my $AggrTag = "";
    $AggrTag = "";
    if ($AggrTag ne "" ){
    
    $mdmDocument->DeselectAll();
    
    for (my $i=0;$i<$nmol;++$i){
        my $mol = $molecules->Item($i);
        my $ID = "";
    #don't seem to be able to access OBJECTID; using a copy for now. This also helps with keeping it as string.
        if($mol->HasProperty("Aggr_OBJECTID")){
        $ID=$mol->GetProperty("Aggr_OBJECTID");
        if( $ID eq $AggrTag) {
            $mol->Select();
            $refMol = $i;
        }
        $mdmDocument->RemoveUserProperty("Aggr_OBJECTID",$mol->ClassName,False);
        }
    
    }
    DSScript::invoke_action("ActionSetSAPReferenceProtein");
    $mdmDocument->DeselectAll();
    }



    $mdmDocument->ShowAggregationSurface(10);


    # Add new Aggr patch groups
    my @AggrRadii = (5,10);
    foreach my $AggrRadius (@AggrRadii) {
    $mdmDocument->FindAggregationSites($AggrRadius);
    }

    #set up data for line plot - per DSC-21034, we want more control over the colouring
    #than could be obtained using 'row index', and FullName doesn't work well if there are multiple molecules.
    my $AggrPlotProperty = "Active Aggr Score";
    my $fullResName = "Residue Specification";
    my $fullChainName = "Chain";
    my $testProperty = "Aggr R" . 10;

    for (my $i=0 ; $i<$nmol ; ++$i){
        my $mol = $molecules->Item($i);
        my $molName = $mol->Name;

        if($molName eq ""){
            $molName = "Mol " . $i;
        }
        my $hasScores = $mol->HasProperty($testProperty);
        if($refMol == -1 && $hasScores){
        $refMol = $i;
        }
        
        my $chains = $mol->AminoAcidChains();
        for (my $j = 0 ;$j < $chains->Count; ++$j){
            my $chain = $chains->Item($j);

            my $chainName = $chain->Name;
            my $fullChain = "";
            if($chainName ne ""){
            $fullChain = $molName . ":" . $chainName;
            $chainName .= ":";
            }else{
                $fullChain = $molName . ":<AminoAcidChain>";
            }

            
    #the commented-out lines set a residue specification which could be used for the line plot -
    #not sure whether we want it so do colouring first
            my $residues = $chain->AminoAcids();
            for (my $k = 0; $k <$residues->Count; ++$k){
                my $residue = $residues->Item($k);
    #			my $resName = $residue->Name;
    #		    my $fullName = $chainName . $resName . ":" . $molName;
                if($hasScores){
    #		      $residue->SetProperty($fullResName,"$fullName");
                $residue->SetProperty($fullChainName, "$fullChain");
                }else{
    #		      $residue->SetProperty($fullResName,"$fullName");
    #			  $residue->SetProperty($fullChainName,"No Aggr Score");
                my $noData = $molName . ": No Scores";
                $residue->SetProperty($fullChainName, "$noData");
                }
            }
        }
    }


    $mdmDocument->UpdateViews();

    # DSScript::invoke_action("ActionShowInTable");



    # my $AggrSiteIDs = "Active Aggr Site ID";

    # #check that patches exist at this radius, otherwise plot can't be created
    #     my $hasPatches = 0;
    #     my $chains = $mdmDocument->AminoAcidChains();
    #     for (my $j = 0 ;$j < $chains->Count; ++$j){
    #         my $chain = $chains->Item($j);
    #         my $residues = $chain->AminoAcids();
    #         for (my $k = 0; $k <$residues->Count; ++$k){
    #             my $residue = $residues->Item($k);
    #             if($residue->HasProperty($AggrSiteIDs)){
    #             $hasPatches = 1;
    #             last;
    #             }
    #         }
    #     }

    # if($hasPatches){
    # my $pointPlot = Chart::Document::CreatePointPlot(
    #     {
    #         MdmDocument      => $mdmDocument,
    #         DataTableTabName => 'AminoAcid',
    #         XAxisProperty    => $AggrSiteIDs,
    #         YAxisProperties  => $AggrPlotProperty
    #     }
    # );
    # # $pointPlot->UpdateViews();
    # }

    # my $linePlot = Chart::Document::CreateLinePlot(
    # {
    #         MdmDocument      => $mdmDocument,
    #         DataTableTabName => 'AminoAcid',
    # #        XAxisProperty    => $fullResName,
    #         XAxisProperty  => "Row Index",
    #         YAxisProperties  => $AggrPlotProperty,
    #         ColorProperty => $fullChainName
    # }
    # );

    # #$linePlot->UpdateViews(); 

    # for(my $i = 0;$i<$nmol; ++$i){
    # my $mol = $molecules->Item($i);
    # if($i == $refMol){
    #     $mol->SetVisible(True);
    # }else{
    #     $mol->SetInvisible(False);
    # }
    # }



#################################### END aggregation score visulization  ###############################################


#07 ################################### 文件IO  ###############################################
    $mdmDocument->UpdateViews();
    $cdrDocument->Save("$work_dir/$basename.bsml",'bsml');
    $mdmDocument->Save("$work_dir/$basename.dsv",'dsv');
    $cdrDocument->Close();
    $sdmDocument->Close();
    $mdmDocument->Close();
    
    rmtree($path_agg_calculation, { verbose => 0, keep_root => 0 });
    rmtree($path_cdrannotation, { verbose => 0, keep_root => 0 });

    my $result1 = join(',',($basename,$totalHPBSAS,$totalY,$totalW,$total_hpb_Y_W));
    my $result2 = join(',',($basename,$v2_totalHPBSAS,$v2_totalY,$v2_totalW,$v2_total_hpb_Y_W));
    return $result1,$result2;
};



open(DATA1,">>$work_dir/HPB_CDR_IMGT_results.csv") || die "HPB_results.txt 文件无法打开, $!";
print DATA1 "ID,HPB_sidechain,Y_sidechain,W_sidechain,total_sidechain\n";
open(DATA2,">>$work_dir/HPB_full_MOL_results.csv") || die "HPB_results.txt 文件无法打开, $!";
print DATA2 "ID,HPB_sidechain,Y_sidechain,W_sidechain,total_sidechain\n";
# 遍历目标目录下的 .pdb 文件
foreach my $file (glob "$work_dir/*pdb") {
    print "Found: $file\n";
    # 在这里处理每个文件
    my ($result1,$result2)= sub_run_calculate_CDR_IMGT_HPB($file,$work_dir);
    # my ($result1,$result2) = ("aaaa","bbbbbb");
    print DATA1 $result1."\n";
    print DATA2 $result2."\n";
};
close(DATA1) || die "无法关闭文件";
close(DATA2) || die "无法关闭文件";