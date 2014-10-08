#' Make sample metadata file for use with IGV.
#'
#' Creates sample metadata file for IGV from a Samplesheet of metadata and FileSheet of file locations.
#'
#'
#' @docType methods
#' @name MakeIGVSampleMetadata
#' @rdname MakeIGVSampleMetadata
#' 
#' @author Thomas Carroll
#'
#' @param SampleSheet A data.frame object with metadata information for samples.
#'    First column must contain unique sample ids. 
#' @param fileSheet A data.frame of file locations. First column must contain the unique sample ids.
#' @param igvdirectory A character of the directory to which sample metadata file is written.
#' @return A character of file location for the IGV sample information file.
#' @import IRanges GenomicRanges XVector Rsamtools tractor.base stringr XML RColorBrewer
#' @include tracktablesFunctions.r
#' @examples
#' 
#' fileLocations <- system.file("extdata",package="tracktables")
#' 
#' bigwigs <- dir(fileLocations,pattern="*.bw",full.names=TRUE)
#' 
#' intervals <- dir(fileLocations,pattern="*.bed",full.names=TRUE)
#' 
#' bigWigMat <- cbind(gsub("_Example.bw","",basename(bigwigs)),
#'                    bigwigs)
#' 
#' intervalsMat <- cbind(gsub("_Peaks.bed","",basename(intervals)),
#'                       intervals)
#' 
#' fileSheet <- merge(bigWigMat,intervalsMat,all=TRUE)
#' 
#' fileSheet <- as.matrix(cbind(fileSheet,NA))
#' 
#' colnames(fileSheet) <- c("SampleName","bigwig","interval","bam")
#' 
#' SampleSheet <- cbind(as.vector(fileSheet[,"SampleName"]),
#'                      c("EBF","H3K4me3","H3K9ac","RNAPol2"),
#'                      c("ProB","ProB","ProB","ProB"))
#' 
#' colnames(SampleSheet) <- c("SampleName","Antibody","Species")
#' MakeIGVSampleMetadata(SampleSheet,fileSheet,igvdirectory=getwd())
#' 
#' @export
MakeIGVSampleMetadata <- function(SampleSheet,fileSheet,igvdirectory){
    write.table("#sampleTable",file.path(igvdirectory,"SampleMetadata.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    colnames(SampleSheet)[1] <- "Linking_id"
    sampleMetadata <- as.matrix(SampleSheet)
    SampleSheet <- as.matrix(fileSheet)
    suppressWarnings(write.table(sampleMetadata,file.path(igvdirectory,"SampleMetadata.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,append=TRUE,sep="\t"))
    BamMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[,"bam"]),"SampleName"],"Bam",sep="_"),SampleSheet[!is.na(SampleSheet[,"bam"]),"SampleName"])
    BigWigMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[,"bigwig"]),"SampleName"],"Bigwig",sep="_"),SampleSheet[!is.na(SampleSheet[,"bigwig"]),"SampleName"])
    IntervalMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[,"interval"]),"SampleName"],"Interval",sep="_"),SampleSheet[!is.na(SampleSheet[,"interval"]),"SampleName"])
    write.table("\n#sampleMapping",file.path(igvdirectory,"SampleMetadata.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,sep="\t")
    write.table("#Bams",file.path(igvdirectory,"SampleMetadata.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,sep="\t")
    write.table(BamMappings,file.path(igvdirectory,"SampleMetadata.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,sep="\t")
    write.table("\n#BigWigs",file.path(igvdirectory,"SampleMetadata.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,sep="\t")
    write.table(BigWigMappings,file.path(igvdirectory,"SampleMetadata.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,sep="\t")
    write.table("\n#Intervals",file.path(igvdirectory,"SampleMetadata.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,sep="\t")
    write.table(IntervalMappings,file.path(igvdirectory,"SampleMetadata.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,sep="\t")
    return(file.path(igvdirectory,"SampleMetadata.txt"))
}

#' Make IGV session XML
#'
#' Creates session XML for IGV from a FileSheet of file locations.
#' 
#'
#'
#' @docType methods
#' @name MakeIGVSessionXML
#' @rdname MakeIGVSessionXML
#' 
#' @author Thomas Carroll
#'
#' @param fileSheet A data.frame of file locations. First column must contain the unique sample ids.
#' @param igvdirectory A character of the directory to which IGV XML session is written.
#' @param XMLname A character of the name for IGV session xml
#' @param genomeName A character of genome for IGV (See IGV user guide for details)
#' @param locusName A character of locus to display in igv on loading (See IGV user guide for details)
#' @param colourBy Character vector of RGB colours to use for colouring displayed BigWigs
#' @param bwScale Character or list of numeric vectors to define scaling type for bigwig. Default is "autoscale".
#' When providing a list, this list must be same length as number of samples and each element have two numeric values corresponding to minimum 
#' and maximum value to be used in setting data range. Currently only "autoscale" or a list of minimum and maximum values are accepted.
#' @return A character of file location for the IGV session XML 
#' @examples
#'  
#' fileLocations <- system.file("extdata",package="tracktables")
#' 
#' bigwigs <- dir(fileLocations,pattern="*.bw",full.names=TRUE)
#' 
#' intervals <- dir(fileLocations,pattern="*.bed",full.names=TRUE)
#' 
#' bigWigMat <- cbind(gsub("_Example.bw","",basename(bigwigs)),
#'                    bigwigs)
#' 
#' intervalsMat <- cbind(gsub("_Peaks.bed","",basename(intervals)),
#'                       intervals)
#' 
#' fileSheet <- merge(bigWigMat,intervalsMat,all=TRUE)
#' 
#' fileSheet <- as.matrix(cbind(fileSheet,NA))
#' 
#' colnames(fileSheet) <- c("SampleName","bigwig","interval","bam")
#' 
#' MakeIGVSessionXML(fileSheet,igvdirectory=getwd(),"Example","mm9")
#' 
#' @export
MakeIGVSessionXML <- function(fileSheet,igvdirectory,XMLname,genomeName,locusName="All",colourBy=NULL,bwScale="autoscale"){
    i <- 1
    SampleSheet <- as.matrix(fileSheet)
    Output <- file.path(igvdirectory,paste(XMLname,".xml",sep=""))
    GlobalNode <- newXMLNode("Global",attrs=c(genome.value=genomeName,groupTracksBy="Linking_id",locus=locusName,version=3))
    ResourcesNode <- newXMLNode("Resources",parent=GlobalNode)
    MetaDataNode <- newXMLNode("Resource",parent=ResourcesNode,attrs=c(name="SampleMetadata",path=relativePath(file.path(igvdirectory,"SampleMetadata.txt"),Output),relativePath=TRUE))
    PanelDataNode <- newXMLNode("Panel",attrs=c(height="350",name="DataPanel",width="1115"),parent=GlobalNode)
    bamFiles <- SampleSheet[,"bam"]
    bigwigFiles <- SampleSheet[,"bigwig"]
    intervalFiles <- SampleSheet[,"interval"]    
    resources <- vector("list")
    if(length(bwScale) == 1){
      bwScale <- rep(bwScale,nrow(SampleSheet))
    }    
    for(i in 1:nrow(SampleSheet)){
        if(!is.null(colourBy)){
            colourIGV <- colourBy[i]
        }else{
            colourIGV <- "0,0,178" 
        }
        if(!is.na(SampleSheet[i,"bam"])){
            NewName <- paste(SampleSheet[i,"SampleName"],"_Bam",sep="")
            resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=relativePath(bamFiles[i],Output),relativePath=TRUE))))
            TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",color="0,0,178",colorOption="UNEXPECTED_PAIR",displayMode="EXPANDED",featureVisibilityWindow="-1",fontSize="10",id=relativePath(bamFiles[i],Output),name=NewName,showDataRange="true",sortByTag="",visible="true"),parent=PanelDataNode)
        }
        if(!is.na(SampleSheet[i,"interval"])){
            NewName <- paste(SampleSheet[i,"SampleName"],"_Interval",sep="")
            resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=relativePath(intervalFiles[i],Output),relativePath=TRUE))))
            TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",color=colourIGV,displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",height="45",id=relativePath(intervalFiles[i],Output),name=NewName,renderer="BASIC_FEATURE",showDataRange="true",sortable="false",visible="true",windowFunction="count"),parent=PanelDataNode)
        }
        if(!is.na(SampleSheet[i,"bigwig"])){
            
            NewName <- paste(SampleSheet[i,"SampleName"],"_Bigwig",sep="")
            resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=relativePath(bigwigFiles[i],Output),relativePath=TRUE))))
            
            if(class(bwScale) == "list"){
                    scaleBigWigIGV <- bwScale[[i]]
                    TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",autoScale="false",color=colourIGV,displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",id=relativePath(bigwigFiles[i],Output),name=NewName,renderer="BAR_CHART",showDataRange="true",visible="true",windowFunction="mean"),parent=PanelDataNode)
                    DisplayRangeNode <-  newXMLNode("DataRange",attrs=c(baseline="0.0",drawBaseline="true",flipAxis="false",maximum=scaleBigWigIGV[2],minimum=scaleBigWigIGV[1],type="LINEAR"),parent=TrackNode)     
            }
            if(bwScale[i] == "autoscale"){
                TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",autoScale="true",color=colourIGV,displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",id=relativePath(bigwigFiles[i],Output),name=NewName,renderer="BAR_CHART",showDataRange="true",visible="true",windowFunction="mean"),parent=PanelDataNode)
                DisplayRangeNode <-  newXMLNode("DataRange",attrs=c(baseline="0.0",drawBaseline="true",flipAxis="false",maximum="50",minimum="5",type="LINEAR"),parent=TrackNode)
                
            }
            if(class(bwScale) != "list" & bwScale[i] != "autoscale"){
                TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",autoScale="false",color=colourIGV,displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",id=relativePath(bigwigFiles[i],Output),name=NewName,renderer="BAR_CHART",showDataRange="true",visible="true",windowFunction="mean"),parent=PanelDataNode)
                DisplayRangeNode <-  newXMLNode("DataRange",attrs=c(baseline="0.0",drawBaseline="true",flipAxis="false",maximum="50",minimum="5",type="LINEAR"),parent=TrackNode)                
            }
            
            
        }
    }  
    saveXML(GlobalNode,file=Output)
  
    return(Output)
}

#' Make HTML pages for IGV sessions (Tracktables Experiment Report)
#'
#' Creates HTML table of sample metadata and all required files for interacting with IGV.
#' 
#'
#'
#' @docType methods
#' @name maketracktable
#' @rdname maketracktable
#' 
#' @author Thomas Carroll
#'
#' @param fileSheet A data frame containing sample file locations (e.g. BigWig locations). 
#' @param SampleSheet A data frame containing sample metadata
#' @param filename Character of name for tracktables HTML report. (.html prefix is added automatically)
#' @param basedirectory Character of directory for tracktables HTML report, IGV sessions and any interval files 
#' @param genome Character of genome for IGV (See IGV user guide for details)
#' @param colourBy Character defining which sample metadata to be used for colouring bigwig files
#' @param bwScale Character or list of numeric vectors to define scaling type for bigwig. Default is "autoscale".
#' When providing a list, this list must be same length as number of samples and each element have two numeric values corresponding to minimum 
#' and maximum value to be used in setting data range. Currently only "autoscale" or a list of minimum and maximum values are accepted.
#' @return An object containing XML document (HTMLInternalDocument,XMLInternalDocument,XMLAbstractDocument) 
#' @examples
#'  
#' fileLocations <- system.file("extdata",package="tracktables")
#' 
#' bigwigs <- dir(fileLocations,pattern="*.bw",full.names=TRUE)
#' 
#' intervals <- dir(fileLocations,pattern="*.bed",full.names=TRUE)
#' 
#' bigWigMat <- cbind(gsub("_Example.bw","",basename(bigwigs)),
#'                    bigwigs)
#' 
#' intervalsMat <- cbind(gsub("_Peaks.bed","",basename(intervals)),
#'                       intervals)
#' 
#' fileSheet <- merge(bigWigMat,intervalsMat,all=TRUE)
#' 
#' fileSheet <- as.matrix(cbind(fileSheet,NA))
#' 
#' colnames(fileSheet) <- c("SampleName","bigwig","interval","bam")
#' 
#' SampleSheet <- cbind(as.vector(fileSheet[,"SampleName"]),
#'                      c("EBF","H3K4me3","H3K9ac","RNAPol2"),
#'                      c("ProB","ProB","ProB","ProB"))
#' 
#' colnames(SampleSheet) <- c("SampleName","Antibody","Species")
#'   HTMLreport <- maketracktable(fileSheet,SampleSheet,
#'                                "IGV_Example.html",
#'                                basedirectory=getwd(),
#'                                "mm9")
#' 
#' @export
maketracktable <- function(fileSheet,SampleSheet,filename,basedirectory,genome,colourBy=NULL,bwScale="autoscale"){
    message("tracktables uses the Datatables javascript libraries.
            For information on Datatables see http://datatables.net/")

    basedirectory <- gsub("/$","",basedirectory)
    MakeIGVSampleMetadata(SampleSheet,fileSheet,basedirectory)
    if(!is.null(colourBy)){
        nOfGroups <- length(unique(SampleSheet[,colourBy]))
        groupColours <- apply(t(col2rgb(brewer.pal(nOfGroups,"Set3"))),1,function(x)paste0(x,collapse=","))[factor(SampleSheet[,colourBy])]
    }else{
        groupColours <- rep("0,0,178",nrow(SampleSheet))
    }
    if(length(bwScale) == 1){
      bwScale <- rep(bwScale,nrow(SampleSheet))
    }
    xmlFiles <- unlist(lapply(seq(1,nrow(fileSheet)),function(x)
        MakeIGVSessionXML(fileSheet[x,,drop=FALSE],
                          basedirectory,
                          paste0(fileSheet[x,1],"igv"),
                          genome,
                          locusName="All",groupColours[x],bwScale[x])
        ))
  
    dataTableJS <- readLines(system.file(package="tracktables","js","datatables.js"))
    jqueryJS <- readLines(system.file(package="tracktables","js","jquery.min.js"))
    dataTableCSS <- readLines(system.file(package="tracktables","js","jquery.datatables.css"))
    dataTableScroller <- readLines(system.file(package="tracktables","js","dataTables.scroller.min.js"))
    tracktablesCSS <- readLines(system.file(package="tracktables","js","tracktables.css"))

    
    giHTMLs <- vector("character",nrow(fileSheet))
    giHTMLLinks <- vector("character",nrow(fileSheet))
    for(l in 1:nrow(fileSheet)){
        if(!is.na(fileSheet[l,"interval"])){
            giHTMLs[l] <- makebedtable(GetGRanges(as.vector(fileSheet[l,"interval"])),paste0(fileSheet[l,"SampleName"],"GI.html"),basedirectory)  
            giHTMLLinks[l] <- paste0("\"<a class=\\\"table\\\" href=\\\"",file.path(basedirectory,basename(giHTMLs[l])),"\\\">Intervals</a>\"")
    }else{
        giHTMLLinks[l] <- shQuote("No Intervals")
      
    }
  }
  
  files <- unlist(lapply(xmlFiles,function(x)relativePath(x,
                                                          gsub("//","/",file.path(basedirectory,filename))
  )))
  t3mp <- "\"<a class=\\\"table\\\" href=\\\"http://localhost:60151/load?file=\".concat(dir.concat(\"/"
  t4mp <- "\\\"\".concat(\""
  t5mp <- "</a>\")))"
  jsMat <- cbind(
    matrix(paste0("\"",as.vector(SampleSheet),"\""),ncol=ncol(SampleSheet),byrow=FALSE),
    paste0(t3mp,files,"&merge=true",t4mp,">",SampleSheet[,1],t5mp),
    giHTMLLinks
  )
  setigv <- paste0("var igvtable = [",paste0(
    "[",apply(jsMat,1,function(x)paste0(
      x,collapse=","))
    ,"]\n",collapse=",")
    ,"];",sep="")
  
  jspart1 <- paste0("var loc = window.location.pathname;\n",
                    "var dir = loc.substring(0, loc.lastIndexOf('/'));\n",setigv,"\n")
  jspart2 <- paste0(
    "$(document).ready(function() {
    $('#demo').html( '<table cellpadding=\"0\" cellspacing=\"0\" border=\"0\" class=\"display\" id=\"example\"></table>' );
    $('#example').dataTable( {
    \"data\": igvtable,\ncolumns:",
    paste0("[",paste0(
      unlist(lapply(c(colnames(SampleSheet),"IGV_Link","Intervals"),function(x)paste0(
        c("{\"title\"",paste0(
          "\"",
          x,"\"}")
        ),collapse=":")
      )),collapse=",\n")
      ,"]")
    ,"\n","} );\n","} );\n")
  jspart1.2 <- paste0(jspart1,jspart2)
  doc <- newXMLDoc(isHTML = TRUE)
  html <- newXMLNode("html",parent=doc)
  head <- newXMLNode("head",parent = html)
  title <- newXMLNode("h2",
                      "Tracktables Report",
                      parent=head)
  css <- newXMLNode("style",
                    attrs=c("style type"="text/css","class"="init"),
                    paste0(dataTableCSS,collapse=""),
                    parent=head)
  tracktablescss <- newXMLNode("style",
                    attrs=c("style type"="text/css","class"="init"),
                    paste0(tracktablesCSS,collapse=""),
                    parent=head)  
  jqueryjs <- newXMLNode("script",
                         attrs=c(type="text/javascript",language="javascript"),
                         paste0(jqueryJS,collapse=""),
                         parent=head)
  datatablejs <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            paste0(dataTableJS,collapse=""),
                            parent=head)
  jspart1.2js <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            jspart1.2,
                            parent=head)
  body <- newXMLNode("body",
                     attrs=c(class="dt-example"),
                     parent=html)
  div <- newXMLNode("div",
                    attrs=c(class="container"),
                    parent=body)
  section <- newXMLNode("section",
                        parent=div)
  divtttext <- newXMLNode("div",
                          attrs=c(id="tttext"),
                          parent=section)
  h1 <- newXMLNode("h1","The Tracktables Experiment Report",
                   parent=divtttext)
  p1 <- newXMLNode("p","
                   This report contains sample information and dynamic links to display and control Broad's Integrative Genome Browser (IGV). This report aims to speed up the organisation and visualisation of genomics data by allowing for the passing of metadata and sample information to IGV and the rapid selection of samples and points of interest using HTML tables.
                   ",
                   parent=divtttext)
  p2 <- newXMLNode("p","Getting started:",
                   parent=divtttext)
  ul1 <- newXMLNode("ul","",
                   parent=divtttext)
  li1 <- newXMLNode("li","To take advantage of the integration with IGV, <b>IGV must be already running </b>on your machine or can be launched now from this <a class=\"main\" href=\"http://www.broadinstitute.org/igv/projects/current/igv.php\">webstart</a>.",
                    parent=ul1,cdata=TRUE)
  li2 <- newXMLNode("li","To load coverage, BAM and/or interval files (bed, narrow peak format etc) simply click the respective sample link in the IGV column.",
                    parent=ul1)
  li3 <- newXMLNode("li","To open a new tracktable containing information on Sample interval files click the link in that sample's repsective Intervals column.",
                    parent=ul1)
  p3 <- newXMLNode("p","For further information on the use of tracktables, please see our Github page or Bioconductor site.",
                   parent=divtttext)
  div2 <- newXMLNode("div",
                     attrs=c(id="demo"),
                     parent=section)
  saveXML(doc,file=file.path(basedirectory,filename),doctype="html")
  return(doc)
}

#' Make HTML pages for interval files or GRanges.
#'
#' Creates HTML pages for interval files or GRanges (Tracktables Interval Report).
#' 
#'
#'
#' @docType methods
#' @name makebedtable
#' @rdname makebedtable
#' 
#' @author Thomas Carroll
#'
#' @param grangesObject A GRanges object.
#' @param name Character of the name for Interval HTML report.
#' @param basedirectory Character of the directory to which HTML report is writen.
#' @return A character of file location for the Tracktables HTML Report
#' @examples
#' data(Intervals)
#' htmlpage <- makebedtable(Intervals,"EBF_PeaksTable.html",getwd())
#' 
#' @export
makebedtable <- function(grangesObject,name,basedirectory){
  
  dataTableJS <- readLines(system.file(package="tracktables","js","datatables.js"))
  jqueryJS <- readLines(system.file(package="tracktables","js","jquery.min.js"))
  dataTableCSS <- readLines(system.file(package="tracktables","js","jquery.datatables.css"))
  dataTableScroller <- readLines(system.file(package="tracktables","js","dataTables.scroller.min.js"))
  tracktablesCSS <- readLines(system.file(package="tracktables","js","tracktables.css"))
  
  grangesFrame <- as.matrix(as.data.frame(grangesObject))
  grangesFrame <- str_trim(grangesFrame)
  jsarray <- paste("[",paste0("[",apply(grangesFrame,1,function(x)paste0(c(shQuote(c(paste0("<a class=\"table\" href=\"http://localhost:60151/goto?locus=",x[1],":",x[2],"-",x[3],"\">IGV</a>"))),shQuote(x)),collapse=",")),"]",collapse=",\n"),"]")
  jsArrayForIGV <- paste0("var igvtable =",jsarray,";\n")
  jspart2 <- paste0(
    "$(document).ready(function() {
    $('#demo').html( '<table cellpadding=\"0\" cellspacing=\"0\" border=\"0\" class=\"display\" id=\"example\"></table>' );
    $('#example').dataTable( {
    deferRender:    true,
    dom:            \"frtiS\",
    scrollY:        200,
    scrollCollapse: true,
    
    \"data\": igvtable,\ncolumns:",
    paste0("[",paste0(
      unlist(lapply(c("IGV_Link",colnames(as.data.frame(grangesObject))),function(x)paste0(
        c("{\"title\"",paste0(
          "\"",
          x,"\"}")
        ),collapse=":")
      )),collapse=",\n")
      ,"]")
    ,"\n","} );\n","} );\n")
  
  jspart1.2 <- paste0(jsArrayForIGV,jspart2)
  doc <- newXMLDoc(isHTML = TRUE)
  html <- newXMLNode("html",parent=doc)
  head <- newXMLNode("head",parent = html)
  title <- newXMLNode("h2",
                      "Tracktables Report",
                      parent=head)
  css <- newXMLNode("style",
                    attrs=c("style type"="text/css","class"="init"),
                    paste0(dataTableCSS,collapse=""),
                    parent=head)
  tracktablescss <- newXMLNode("style",
                               attrs=c("style type"="text/css","class"="init"),
                               paste0(tracktablesCSS,collapse=""),
                               parent=head)  
  jqueryjs <- newXMLNode("script",
                         attrs=c(type="text/javascript",language="javascript"),
                         paste0(jqueryJS,collapse=""),
                         parent=head)
  datatablejs <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            paste0(dataTableJS,collapse=""),
                            parent=head)
  datatableScroller <- newXMLNode("script",
                                  attrs=c(type="text/javascript",language="javascript"),
                                  paste0(dataTableScroller,collapse=""),
                                  parent=head)  
  jspart1.2js <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            jspart1.2,
                            parent=head)
  body <- newXMLNode("body",
                     attrs=c(class="dt-example"),
                     parent=html)
  div <- newXMLNode("div",
                    attrs=c(class="container"),
                    parent=body)
  section <- newXMLNode("section",
                        parent=div)
  divtttext <- newXMLNode("div",
                     attrs=c(id="tttext"),
                     parent=section)
  h1 <- newXMLNode("h1","The Tracktables Interval Report",
                   parent=divtttext)
  p1 <- newXMLNode("p","
                   This report contains genomic interval coordinates,  metadata and dynamic links to control the region displayed within Broad's Integrative Genome Browser (IGV). This alows fort rapid visualisation and interrogation of points of interest within the Genome Browser using HTML tables.
                   ",
                   parent=divtttext)
  p2 <- newXMLNode("p","Getting started:",
                   parent=divtttext)
  ul1 <- newXMLNode("ul","",
                    parent=divtttext)
  li1 <- newXMLNode("li","To take advantage of the integration with IGV, <b>IGV must be already running </b>on your machine or can be launched now from this <a class=\"main\" href=\"http://www.broadinstitute.org/igv/projects/current/igv.php\">webstart</a>.",
                    parent=ul1,cdata=TRUE)
  li2 <- newXMLNode("li","To change IGV display to the region of interest, simply click the respective Interval link in the IGV column.",
                    parent=ul1)
  p3 <- newXMLNode("p","For further information on the use of tracktables, please see our Github page or Bioconductor site.",
                   parent=divtttext)  
  div2 <- newXMLNode("div",
                     attrs=c(id="demo"),
                     parent=section)
  saveXML(doc,file=file.path(basedirectory,name),doctype="html")
  
  
}



GetGRanges <- function(LoadFile,AllChr=NULL,ChrOfInterest=NULL,simple=FALSE,sepr="\t",simplify=FALSE){
  
  if(class(LoadFile) == "GRanges"){
    RegionRanges <- LoadFile
    if(simplify){
      RegionRanges <- GRanges(seqnames(RegionRanges),ranges(RegionRanges))
    }
  }else{
    if(class(LoadFile) == "character"){
      RangesTable <- read.delim(LoadFile,sep=sepr,header=FALSE,comment.char="#")
    }else if(class(LoadFile) == "matrix"){
      RangesTable <- as.data.frame(LoadFile)
    } else{
      RangesTable <- as.data.frame(LoadFile)
    }
    Chromosomes <- as.vector(RangesTable[,1])
    Start <- as.numeric(as.vector(RangesTable[,2]))
    End <- as.numeric(as.vector(RangesTable[,3]))
    RegionRanges <- GRanges(seqnames=Chromosomes,ranges=IRanges(start=Start,end=End))
    if(simple == FALSE){
      if(ncol(RangesTable) > 4){
        ID <- as.vector(RangesTable[,4])
        Score <- as.vector(RangesTable[,5])
        if(ncol(RangesTable) > 6){
          Strand <- rep("*",nrow(RangesTable))
          RemainderColumn <- as.data.frame(RangesTable[,-c(1:6)])
          elementMetadata(RegionRanges) <- cbind(ID,Score,Strand,RemainderColumn)
        }else{
          elementMetadata(RegionRanges) <- cbind(ID,Score)
        }
      }
    }
  } 
  return(RegionRanges)
}

#' Make IGV session XML and sample information file
#'
#' #' Creates IGV session XML and sample information file from a Samplesheet of metadata and FileSheet of file locations.
#' 
#'
#'
#' @docType methods
#' @name MakeIGVSession
#' @rdname MakeIGVSession
#' 
#' @author Thomas Carroll
#'
#' @param SampleSheet A data.frame object with metadata information for samples.
#'  First column must contain unique sample ids. 
#' @param fileSheet A data.frame of file locations. First column must contain the unique sample ids.
#' @param igvdirectory A character of the directory to which sample metadata file is written.
#' @param XMLname A character of the name for IGV session xml
#' @param genomeName A character of genome for IGV (See IGV user guide for details)
#' @param locusName A character of locus to display in igv on loading (See IGV user guide for details)
#' @param colourBy Character defining which sample metadata to be used for colouring bigwig files
#' @param bwScale Character or list of numeric vectors to define scaling type for bigwig. Default is "autoscale".
#' When providing a list, this list must be same length as number of samples and each element have two numeric values corresponding to minimum 
#' and maximum value to be used in setting data range. Currently only "autoscale" or a list of minimum and maximum values are accepted.
#' @return A character of file location for the IGV session XML
#' @examples
#'  
#' fileLocations <- system.file("extdata",package="tracktables")
#' 
#' bigwigs <- dir(fileLocations,pattern="*.bw",full.names=TRUE)
#' 
#' intervals <- dir(fileLocations,pattern="*.bed",full.names=TRUE)
#' 
#' bigWigMat <- cbind(gsub("_Example.bw","",basename(bigwigs)),
#'                    bigwigs)
#' 
#' intervalsMat <- cbind(gsub("_Peaks.bed","",basename(intervals)),
#'                       intervals)
#' 
#' fileSheet <- merge(bigWigMat,intervalsMat,all=TRUE)
#' 
#' fileSheet <- as.matrix(cbind(fileSheet,NA))
#' 
#' colnames(fileSheet) <- c("SampleName","bigwig","interval","bam")
#' 
#' SampleSheet <- cbind(as.vector(fileSheet[,"SampleName"]),
#'                      c("EBF","H3K4me3","H3K9ac","RNAPol2"),
#'                      c("ProB","ProB","ProB","ProB"))
#' 
#' colnames(SampleSheet) <- c("SampleName","Antibody","Species")
#' MakeIGVSession(SampleSheet,fileSheet,igvdirectory=getwd(),"Example","mm9")
#' 
#' @export
MakeIGVSession <- function(SampleSheet,fileSheet,igvdirectory,XMLname,genomeName,locusName="All",colourBy=NULL,bwScale="autoscale"){
  if(!is.null(colourBy)){
    nOfGroups <- length(unique(SampleSheet[,colourBy]))
    groupColours <- apply(t(col2rgb(brewer.pal(nOfGroups,"Set3"))),1,function(x)paste0(x,collapse=","))[factor(SampleSheet[,colourBy])]
  }else{
    groupColours <- rep("0,0,178",nrow(SampleSheet))
  }
  
  MakeIGVSampleMetadata(SampleSheet,fileSheet,igvdirectory)
  sessionxml <- MakeIGVSessionXML(fileSheet,igvdirectory,XMLname,genomeName,locusName="All",colourBy=groupColours,bwScale=bwScale)  
  return(sessionxml)
}

#' Example genomic intervals
#'
#' This dataset contains peaks from an in-house EBF1 ChIP-seq 
#'
#' \itemize{
#' \item Intervals GRanges object containing EBF1 peaks
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Intervals
#' @usage data(Intervals)
#' @return A GRanges object with two rows
NULL
