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
#' @import IRanges GenomicRanges XVector Rsamtools tractor.base stringr XML RColorBrewer methods
#' @include tracktablesFunctions.R
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
#' @param igvParams An object of class igvParam containing display parameters for IGV.
#' When providing a list, this list must be same length as number of samples and each element have two numeric values corresponding to minimum 
#' and maximum value to be used in setting data range. Currently only "autoscale" or a list of minimum and maximum values are accepted.
#' @param writedirectory A character of the directory to which files will be written. Default is set
#' to igvdirectory argument. 
#' @param full.xml.paths Boolean of whether reference to XML and sample information files should 
#' be by relative or absolute paths. Default is FALSE
#' @param full.file.paths Boolean of whether reference to sample bigWig/Bam/interval files should 
#' be by relative or absolute paths. Default is FALSE
#' @param use.path.asis Boolean of whether paths to files in samplesheet should be used as is. 
#' If TRUE overrides sample files paths specified by full.file.paths argument. Default is FALSE
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
MakeIGVSessionXML <- function(fileSheet,igvdirectory,XMLname,genomeName,locusName="All",
                              colourBy=NULL,igvParams=igvParam(),writedirectory=NULL,
                              full.xml.paths=FALSE,full.file.paths=FALSE,use.path.asis=FALSE){
    i <- 1
    if(is.null(writedirectory)){
      writedirectory <- igvdirectory
    }
    SampleSheet <- as.matrix(fileSheet)
    if(class(igvParams) == "igvParam"){
      igvParams <- rep(list(igvParams),nrow(fileSheet))
    }
    if(!full.xml.paths){
      Output <- file.path(writedirectory,paste(XMLname,".xml",sep=""))
    }else{
      Output <- file.path(writedirectory,paste(XMLname,".xml",sep=""))
    }
    GlobalNode <- newXMLNode("Global",attrs=c(genome.value=genomeName,groupTracksBy="Linking_id",locus=locusName,version=3))
    ResourcesNode <- newXMLNode("Resources",parent=GlobalNode)
    if(full.xml.paths){
      sampleMetadataPath <- file.path(igvdirectory,"SampleMetadata.txt")
      relativePathSampleMetadataFlag <- "FALSE"
    }else{      
      sampleMetadataPath <- relativePath(file.path(igvdirectory,"SampleMetadata.txt"),Output)
      relativePathSampleMetadataFlag <- "TRUE"
      if(full.file.paths){
        sampleMetadataPath <- relativePath(file.path(writedirectory,"SampleMetadata.txt"),Output)  
      }
    }
    MetaDataNode <- newXMLNode("Resource",parent=ResourcesNode,attrs=c(name="SampleMetadata",path=sampleMetadataPath,relativePath=relativePathSampleMetadataFlag))
    PanelDataNode <- newXMLNode("Panel",attrs=c(height="350",name="DataPanel",width="1115"),parent=GlobalNode)
    bamFiles <- SampleSheet[,"bam"]
    bigwigFiles <- SampleSheet[,"bigwig"]
    intervalFiles <- SampleSheet[,"interval"]    
    resources <- vector("list")
    for(i in 1:nrow(SampleSheet)){
        if(!is.null(colourBy)){
            colourIGVbam <- colourBy[i]
            colourIGVbigWig <- colourBy[i]
            colourIGVinterval <- colourBy[i]
        }else{
            #print(paste0(col2rgb(igvParams[[i]]@bigwig.color),collapse=","))
            colourIGVbam <- paste0(col2rgb(igvParams[[i]]@bam.color),collapse=",")
            colourIGVbigWig <- paste0(col2rgb(igvParams[[i]]@bigwig.color),collapse=",")
            colourIGVinterval <- paste0(col2rgb(igvParams[[i]]@interval.color),collapse=",")              
        }
        if(!is.na(SampleSheet[i,"bam"])){
            NewName <- paste(SampleSheet[i,"SampleName"],"_Bam",sep="")
            if(full.file.paths){
              bamFilePath <-  file.path(igvdirectory,basename(unname(bamFiles[i])))
              relativePathBamFlag <- FALSE
            }else{
              bamFilePath <- relativePath(bamFiles[i],Output)
              relativePathBamFlag <- TRUE
            }
            if(use.path.asis){
              bamFilePath <- unname(bamFiles[i])
              relativePathBamFlag <- FALSE              
            }            
            resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=bamFilePath,relativePath=relativePathBamFlag))))
            TrackNode <-  newXMLNode("Track",attrs=c(id=bamFilePath,
                                                     name=NewName,
                                                     showDataRange="true",
                                                     color=colourIGVbam,
                                                     altColor=paste0(col2rgb(igvParams[[i]]@bam.altColor),collapse=","),
                                                     colorOption=igvParams[[i]]@bam.colorOption,
                                                     displayMode=igvParams[[i]]@bam.displayMode,                                                     
                                                     featureVisibilityWindow=as.character(igvParams[[i]]@bam.featureVisibilityWindow),
                                                     fontSize=as.character(igvParams[[i]]@bam.fontSize),                                                     
                                                     autoScale=igvParams[[i]]@bam.autoScale,
                                                     showSpliceJunctions=igvParams[[i]]@bam.showSpliceJunctions,
                                                     flagUnmappedPairs=igvParams[[i]]@bam.flagUnmappedPairs,
                                                     colorByTag=igvParams[[i]]@bam.colorByTag,
                                                     groupByTag=igvParams[[i]]@bam.groupByTag,
                                                     sortByTag=igvParams[[i]]@bam.sortByTag,
                                                     minInsertSize=igvParams[[i]]@bam.minInsertSize,
                                                     maxInsertSize=igvParams[[i]]@bam.maxInsertSize,
                                                     shadeBasesOption=igvParams[[i]]@bam.shadeBasesOption,
                                                     shadeCenters=igvParams[[i]]@bam.shadeCenters,
                                                     showAllBases=igvParams[[i]]@bam.showAllBases,
                                                     visible="true"),
                                     parent=PanelDataNode)
        }
        if(!is.na(SampleSheet[i,"interval"])){
            NewName <- paste(SampleSheet[i,"SampleName"],"_Interval",sep="")
            if(full.file.paths){
              intervalFilePath <- file.path(igvdirectory,basename(unname(intervalFiles[i])))
              relativePathIntervalFlag <- FALSE
            }else{
              intervalFilePath <- relativePath(intervalFiles[i],Output)
              relativePathIntervalFlag <- TRUE
            }
            if(use.path.asis){
              intervalFilePath <- unname(intervalFiles[i])
              relativePathIntervalFlag <- FALSE              
            }
            resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=intervalFilePath,relativePath=relativePathIntervalFlag))))
            TrackNode <-  newXMLNode("Track",
                                     attrs=c(id=intervalFilePath,
                                             name=NewName,
                                             altColor=paste0(col2rgb(igvParams[[i]]@interval.altColor),collapse=","),
                                             color=colourIGVinterval,
                                             displayMode=igvParams[[i]]@interval.displayMode,
                                             featureVisibilityWindow=as.character(igvParams[[i]]@interval.featureVisibilityWindow),
                                             fontSize=as.character(igvParams[[i]]@interval.fontSize),
                                             height=as.character(igvParams[[i]]@interval.height),
                                             renderer=igvParams[[i]]@interval.renderer,
                                             showDataRange="true",
                                             sortable=igvParams[[i]]@interval.sortable,
                                             visible=igvParams[[i]]@interval.visible,
                                             windowFunction=igvParams[[i]]@interval.windowFunction,
                                             autoScale=igvParams[[i]]@interval.autoScale,
                                             normalize=igvParams[[i]]@interval.normalize                                             
                                             ),
                                     parent=PanelDataNode)
        }
        if(!is.na(SampleSheet[i,"bigwig"])){
            
            NewName <- paste(SampleSheet[i,"SampleName"],"_Bigwig",sep="")
            if(full.file.paths){
              bigWigFilePath <- file.path(igvdirectory,basename(unname(bigwigFiles[i])))
              relativePathbigWigFlag <- FALSE
            }else{
              bigWigFilePath <- relativePath(bigwigFiles[i],Output)
              relativePathbigWigFlag <- TRUE
            }
            if(use.path.asis){
              bigWigFilePath <- unname(bigwigFiles[i])
              relativePathbigWigFlag <- FALSE              
            }             
            resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=bigWigFilePath,relativePath=relativePathbigWigFlag))))
            TrackNode <-  newXMLNode("Track",attrs=c(id=bigWigFilePath,
                                                     name=NewName,
                                                     autoScale=igvParams[[i]]@bigwig.autoScale,
                                                     altColor=paste0(col2rgb(igvParams[[i]]@bigwig.altColor),collapse=","),
                                                     color=colourIGVbigWig,
                                                     displayMode=igvParams[[i]]@bigwig.displayMode,
                                                     featureVisibilityWindow=as.character(igvParams[[i]]@bigwig.featureVisibilityWindow),
                                                     fontSize=as.character(igvParams[[i]]@bigwig.fontSize),
                                                     renderer=igvParams[[i]]@bigwig.renderer,
                                                     showDataRange="true",
                                                     visible=igvParams[[i]]@bigwig.visible,
                                                     sortable=igvParams[[i]]@bigwig.sortable,
                                                     normalize=igvParams[[i]]@bigwig.normalize,
                                                     windowFunction=igvParams[[i]]@bigwig.windowFunction),
                                     parent=PanelDataNode)
            DisplayRangeNode <-  newXMLNode("DataRange",attrs=c(baseline=igvParams[[i]]@bigwig.baseline,
                                                                drawBaseline=igvParams[[i]]@bigwig.drawBaseline,
                                                                flipAxis=igvParams[[i]]@bigwig.flipAxis,
                                                                maximum=igvParams[[i]]@bigwig.maximum,
                                                                minimum=igvParams[[i]]@bigwig.minimum,
                                                                type=igvParams[[i]]@bigwig.type),
                                            parent=TrackNode)                           
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
#' @param igvParams An object of class igvParam containing display parameters for IGV.
#' When providing a list, this list must be same length as number of samples and each element have two numeric values corresponding to minimum 
#' and maximum value to be used in setting data range. Currently only "autoscale" or a list of minimum and maximum values are accepted.
#' @param writedirectory A character of the directory to which files will be written. Default is set
#' to igvdirectory argument. 
#' @param full.xml.paths Boolean of whether reference to XML and sample information files should 
#' be by relative or absolute paths. Default is FALSE
#' @param full.file.paths Boolean of whether reference to sample bigWig/Bam/interval files should 
#' be by relative or absolute paths. Default is FALSE
#' @param use.path.asis Boolean of whether paths to files in samplesheet should be used as is. 
#' If TRUE overrides sample files paths specified by full.file.paths argument. Default is FALSE
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
maketracktable <- function(fileSheet,SampleSheet,filename,basedirectory,genome,colourBy=NULL,
                           igvParams=igvParam(),writedirectory=NULL,
                           full.xml.paths=FALSE,full.file.paths=FALSE,use.path.asis=FALSE){
      message("tracktables uses the Datatables javascript libraries.
            For information on Datatables see http://datatables.net/")
    if(class(igvParams) == "igvParam"){
      igvParams <- rep(list(igvParams),nrow(fileSheet))
    }
    if(class(igvParams) == "list"){
      if(length(igvParams) != nrow(fileSheet)){
        igvParams <- igvParams[1]
        message("igvParams is not the same length as fileSheet. Only first igvParams in list will be used")
      }
    }
    
    basedirectory <- gsub("/$","",basedirectory)
    if(is.null(writedirectory)) writedirectory <- basedirectory
    
    if(!full.xml.paths){
      MakeIGVSampleMetadata(SampleSheet,fileSheet,writedirectory)
    }else{
      MakeIGVSampleMetadata(SampleSheet,fileSheet,writedirectory)
    }
    if(!is.null(colourBy)){
        nOfGroups <- length(unique(SampleSheet[,colourBy]))
        groupColours <- apply(t(col2rgb(brewer.pal(nOfGroups,"Set3"))),1,function(x)paste0(x,collapse=","))[factor(SampleSheet[,colourBy])]
    }else{
        groupColours <- NULL
    }
    xmlFiles <- unlist(lapply(seq(1,nrow(fileSheet)),function(x)
        MakeIGVSessionXML(fileSheet[x,,drop=FALSE],
                          basedirectory,
                          paste0(fileSheet[x,1],"igv"),
                          genome,
                          locusName="All",groupColours[x],igvParams[x],
                          writedirectory=writedirectory,
                          full.xml.paths=full.xml.paths,full.file.paths=full.file.paths,
                          use.path.asis=use.path.asis)
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
            if(!full.xml.paths){
            giHTMLs[l] <- makebedtable(GetGRanges(as.vector(fileSheet[l,"interval"])),paste0(fileSheet[l,"SampleName"],"GI.html"),writedirectory)  
            htmlfile <- unlist(lapply(giHTMLs[l],function(x)relativePath(x, gsub("//","/",file.path(writedirectory,filename)))))
            giHTMLLinks[l] <- paste0("\"<a class=\\\"table\\\" href=\\\"",
                                     htmlfile,
                                     "\\\">Intervals</a>\"")
            }else{
              giHTMLs[l] <- makebedtable(GetGRanges(as.vector(fileSheet[l,"interval"])),paste0(fileSheet[l,"SampleName"],"GI.html"),writedirectory)  
              giHTMLLinks[l] <- paste0("\"<a class=\\\"table\\\" href=\\\"",file.path(basedirectory,basename(giHTMLs[l])),"\\\">Intervals</a>\"") 
            }
    }else{
        giHTMLLinks[l] <- shQuote("No Intervals")
      
    }
  }

  if(!full.xml.paths){
    if(full.file.paths){
      message(paste0("full.file.paths is set to true,\nrelative paths will be created to XML from writedirectory: ",writedirectory))
    files <- unlist(lapply(xmlFiles,function(x)relativePath(x,
                                                            gsub("//","/",file.path(writedirectory,filename)))))
                                                              
    }else{
    files <- unlist(lapply(xmlFiles,function(x)relativePath(x,
                                                            gsub("//","/",file.path(basedirectory,filename)))))
    }                                                        
    t3mp <- "\"<a class=\\\"table\\\" href=\\\"http://localhost:60151/load?file=\".concat(dir.concat(\"/"
    t4mp <- "\\\"\".concat(\""
    t5mp <- "</a>\")))"
    jsMat <- cbind(
      matrix(paste0("\"",as.vector(SampleSheet),"\""),ncol=ncol(SampleSheet),byrow=FALSE),
      paste0(t3mp,files,"&merge=true",t4mp,">",SampleSheet[,1],t5mp),
      giHTMLLinks
    )
  }else{
    files <- unlist(xmlFiles)
    files <- file.path(basedirectory,basename(files))
    t3mp <- "\"<a class=\\\"table\\\" href=\\\"http://localhost:60151/load?file="
    t4mp <- ""
    t5mp <- "</a>\""
    jsMat <- cbind(
      matrix(paste0("\"",as.vector(SampleSheet),"\""),ncol=ncol(SampleSheet),byrow=FALSE),
      paste0(t3mp,files,"&merge=true\\\"",t4mp,">",SampleSheet[,1],t5mp),
      giHTMLLinks
    )
    
  }
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
  if(!full.xml.paths){
    saveXML(doc,file=file.path(writedirectory,filename),doctype="html")
  }else{
    saveXML(doc,file=file.path(writedirectory,filename),doctype="html")
  }
  newdoc <- newXMLDoc(isHTML = T)
  html <- newXMLNode("html",parent=newdoc)
  div <- newXMLNode("head",
                    parent=html)
  #css <- newXMLNode("style",
  #                  attrs=c("style type"="text/css","class"="init"),
  #                          paste0(dataTableCSS,collapse=""),
  #                  parent=div)
  #tracktablescss <- newXMLNode("style",
  #                             attrs=c("style type"="text/css","class"="init"),
  #                                     paste0(tracktablesCSS,collapse=""),
  #                             parent=div)  
  #jqueryjs <- newXMLNode("script",
  #                       attrs=c(type="text/javascript",language="javascript"),
  #                               paste0(jqueryJS,collapse=""),
  #                       parent=div)
  #datatablejs <- newXMLNode("script",
  #                          attrs=c(type="text/javascript",language="javascript"),
  #                                  paste0(dataTableJS,collapse=""),
  #                          parent=div)
  #jspart1.2js <- newXMLNode("script",
  #                          attrs=c(type="text/javascript",language="javascript"),
  #                          paste0(jspart1.2,collapse=""),
  #                          parent=div)
  #div2 <- newXMLNode("body",
  #                  parent=html)
  #section <- newXMLNode("section",
  #                      parent=div2)
  #divtttext <- newXMLNode("div",
  #                        attrs=c(id="tttext"),
  #                        parent=section)

  #ul1 <- newXMLNode("ul","",
  #                  parent=divtttext)
  #li1 <- newXMLNode("li","To take advantage of the integration with IGV, <b>IGV must be already running </b>on your machine or can be launched now from this <a class=\"main\" href=\"http://www.broadinstitute.org/igv/projects/current/igv.php\">webstart</a>.",
  #                  parent=ul1,cdata=TRUE)
  #div2 <- newXMLNode("div",
  #                   attrs=c(id="demo"),
  #                   parent=section)

  css <- newXMLNode("link",
                    attrs=c("rel"="stylesheet","type"="text/css",
                            "href"="https://cdn.rawgit.com/ThomasCarroll/tracktables-Data/master/js/jquery.datatables.css"),
                    parent=div)
  tracktablescss <- newXMLNode("link",
                               attrs=c("rel"="stylesheet","type"="text/css",
                               "href"="https://cdn.rawgit.com/ThomasCarroll/tracktables-Data/master/js/tracktables.css"),
                               parent=div)  
  jqueryjs <- newXMLNode("script",
                         attrs=c(type="text/javascript",language="javascript",
                                "src"="https://cdn.rawgit.com/ThomasCarroll/tracktables-Data/master/js/jquery.min.js"),
                         parent=div)
  datatablejs <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript",
                            "src"="https://cdn.rawgit.com/ThomasCarroll/tracktables-Data/master/js/datatables.js"),
                            parent=div)
  jspart1.2js <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            paste0(jspart1.2,collapse=""),
                            parent=div)
  div2 <- newXMLNode("body",
                    parent=html)
  section <- newXMLNode("section",
                        parent=div2)
  divtttext <- newXMLNode("div",
                          attrs=c(id="tttext"),
                          parent=section)
  
  ul1 <- newXMLNode("ul","",
                    parent=divtttext)
  li1 <- newXMLNode("li","To take advantage of the integration with IGV, <b>IGV must be already running </b>on your machine or can be launched now from this <a class=\"main\" href=\"http://www.broadinstitute.org/igv/projects/current/igv.php\">webstart</a>.",
                    parent=ul1,cdata=TRUE)
  div2 <- newXMLNode("div",
                     attrs=c(id="demo"),
                     parent=section)
  
  
  return(
    gsub("\\\\\"","'",
      gsub("</body>|<body>|</html>|<html>|</head>|<head>","",
              saveXML(newdoc))
    )
  )
  
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
  grangesFrame <- apply(grangesFrame,2,str_trim)
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
#' @param igvParams An object of class igvParam containing display parameters for IGV.
#' When providing a list, this list must be same length as number of samples and each element have two numeric values corresponding to minimum 
#' and maximum value to be used in setting data range. Currently only "autoscale" or a list of minimum and maximum values are accepted.
#' @param writedirectory A character of the directory to which files will be written. Default is set
#' to igvdirectory argument. 
#' @param full.xml.paths Boolean of whether reference to XML and sample information files should 
#' be by relative or absolute paths. Default is FALSE
#' @param full.file.paths Boolean of whether reference to sample bigWig/Bam/interval files should 
#' be by relative or absolute paths. Default is FALSE
#' @param use.path.asis Boolean of whether paths to files in samplesheet should be used as is. 
#' If TRUE overrides sample files paths specified by full.file.paths argument. Default is FALSE
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
MakeIGVSession <- function(SampleSheet,fileSheet,igvdirectory,XMLname,genomeName,locusName="All",colourBy=NULL,
                           igvParams=igvParam(),writedirectory=NULL,full.xml.paths=FALSE,full.file.paths=FALSE,
                           use.path.asis=FALSE){
  if(!is.null(colourBy)){
    nOfGroups <- length(unique(SampleSheet[,colourBy]))
    groupColours <- apply(t(col2rgb(brewer.pal(nOfGroups,"Set3"))),1,function(x)paste0(x,collapse=","))[factor(SampleSheet[,colourBy])]
  }else{
    groupColours <- NULL
  }
  
  MakeIGVSampleMetadata(SampleSheet,fileSheet,igvdirectory)
  sessionxml <- MakeIGVSessionXML(fileSheet,igvdirectory,XMLname,genomeName,locusName="All",colourBy=groupColours,igvParams=igvParams,full.xml.paths=full.xml.paths,full.file.paths=full.file.paths,use.path.asis=use.path.asis)  
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

#' Parameters for displaying bigwigs, bams and intervals in IGV
#'
#' Use igvParam( object)) to create a parameter object to control IGV display invoked from maketracktable() report or 
#' from makeIGVSession() XML files. Parameters for bigwig, bam and intervals files may be provided. See IGV manual for a full
#' description of options.
#'
#' @docType class
#' @rdname igvParam
#' @aliases igvParam-class igvParam
#' @author Thomas Carroll
#' @param bigwig.altColor A character vector of alternate colour for bigwigs displayed in IGV.
#' @param bigwig.color A character vector of main colour for bigwigs displayed in IGV.
#' @param bigwig.displayMode A character vector of display mode for bigwigs displayed in IGV.
#' @param bigwig.featureVisibilityWindow A numeric vector of feature visibility window for bigwigs displayed in IGV (Defaut -1).
#' @param bigwig.fontSize A numeric vector of font size for bigwigs displayed in IGV.
#' @param bigwig.autoScale A character vector ("true"/"false") to indicate whether data is autoscaled for bigwigs displayed in IGV (Default "true").
#' @param bigwig.normalize A character vector ("true"/"false") to indicate whether data is normalised for bigwigs displayed in IGV (Default "false").
#' @param bigwig.renderer A character vector of renderer for bigwigs displayed in IGV (Default "BAR_CHART").
#' @param bigwig.sortable A character vector ("true"/"false") to indicate whether data is sortable for bigwigs displayed in IGV (Default "true").
#' @param bigwig.visible A character vector ("true"/"false") to indicate whether data is visible for bigwigs displayed in IGV (Default "true").
#' @param bigwig.windowFunction A character vector ("true"/"false") of window function for bigwigs displayed in IGV (Default "mean").
#' @param bigwig.baseline  A numeric vector of baseline bigwigs displayed in IGV.
#' @param bigwig.drawBaseline  A character vector ("true"/"false") of whether to draw baseline for bigwigs displayed in IGV (Default "true").
#' @param bigwig.flipAxis  A character vector ("true"/"false") to indicate whether to flip axis for bigwigs displayed in IGV (Default "false").
#' @param bigwig.maximum  A numeric vector of maximum value to display for bigwigs displayed in IGV (Default 50).
#' @param bigwig.minimum A numeric vector of minimum value to display for bigwigs displayed in IGV (Default 0).
#' @param bigwig.type A character vector of display type for bigwigs displayed in IGV (Default "LINEAR").
#' @param interval.altColor  A character vector of alternate colour for intervals displayed in IGV. 
#' @param interval.color  A character vector of main colour for intervals displayed in IGV.
#' @param interval.autoScale A character vector ("true"/"false") to indicate whether data is autoscaled for intervals displayed in IGV (Default "false").
#' @param interval.displayMode A character vector of display mode for intervals displayed in IGV (Default "COLLAPSED").                            
#' @param interval.featureVisibilityWindow A numeric vector of feature visibility window for intervals displayed in IGV (Defaut -1).
#' @param interval.fontSize  A numeric vector of font size for intervals displayed in IGV.
#' @param interval.height  A numeric vector of height for intervals displayed in IGV.
#' @param interval.normalize A character vector ("true"/"false") to indicate whether data is normalised for intervals displayed in IGV (Default "false").
#' @param interval.renderer A character vector of renderer for intervals displayed in IGV (Default "BASIC_FEATURE").
#' @param interval.sortable A character vector ("true"/"false") to indicate whether data is sortable for intervals displayed in IGV (Default "true").
#' @param interval.visible A character vector ("true"/"false") to indicate whether data is visible for intervals displayed in IGV (Default "true").
#' @param interval.windowFunction  A character vector ("true"/"false") of window function for intervals displayed in IGV (Default "count").
#' @param bam.altColor A character vector of alternate colour for bam files displayed in IGV. 
#' @param bam.color  A character vector of main colour for bam files displayed in IGV. 
#' @param bam.autoScale A character vector ("true"/"false") to indicate whether data is autoscaled for bam files displayed in IGV (Default "false").
#' @param bam.displayMode A character vector of display mode for bam files displayed in IGV (Default "EXPANDED").                           
#' @param bam.featureVisibilityWindow  A numeric vector of feature visibility window for bam files displayed in IGV (Defaut -1).
#' @param bam.fontSize A numeric vector of font size for intervals displayed in IGV.
#' @param bam.showSpliceJunctions A character vector ("true"/"false") to indicate whether to show splice juntions for bam files displayed in IGV (Default "false").
#' @param bam.colorByTag  A character vector to indicate whether to colour reads by Tags for Bam files (Defualt = "").
#' @param bam.colorOption A character vector of option to highlight Tags for Bam files (Defualt = "UNEXPECTED_PAIR").
#' @param bam.flagUnmappedPairs A character vector ("true"/"false") to indicate whether to flaf unmapped pairs for bam files displayed in IGV (Default "false").
#' @param bam.groupByTag A character vector ("true"/"false") to indicate how to groups reads by Tag for bam files displayed in IGV (Default "").
#' @param bam.maxInsertSize A numeric vector of maximum insert size to display for Bam files.
#' @param bam.minInsertSize A numeric vector of minimum insert size to display for Bam files.
#' @param bam.shadeBasesOption A character vector of option to shade bases for Bam files (Default "QUALITY").
#' @param bam.shadeCenters A character vector ("true"/"false") to indicate whether to shade centres for bam files displayed in IGV (Default "false").
#' @param bam.showAllBases A character vector ("true"/"false") to indicate to show all bases for bam files displayed in IGV (Default "false").
#' @param bam.sortByTag A character vector ("true"/"false") to indicate how to sort reads by Tag for bam files displayed in IGV (Default "").
#' @return An igvParam class object to use with maketracktable,MakeIGVSession and makeIGVSessionMXL
#' @examples
#' 
#' ## Simple initialisation of an IGVParam object  
#' igvDisplayParams <- igvParam()
#' 
#' 
#' ## More custom initialisation of an IGVParam object .
#' igvDisplayParams <- igvParam(bigwig.color="red",bigwig.autoScale = "false",
#' bigwig.minimum = 10,bigwig.maximum = 100)
#' 
#' 
#' # See full parameters and IGV online manual for more details on customistions  
#'  
#'  
#' ## Use igvParams with maketracktables function to customise bigwig display colour and data range.
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
#' MakeIGVSession(SampleSheet,fileSheet,
#' igvdirectory=getwd(),"Example","mm9",
#' igvParams=igvDisplayParams)
#' 
igvParam <- setClass("igvParam",
                     slots = c(bigwig.altColor="character", 
                               bigwig.color="character",
                               bigwig.autoScale="character",
                               bigwig.displayMode="character",
                               bigwig.featureVisibilityWindow="numeric",
                               bigwig.fontSize="numeric",
                               bigwig.normalize="character",
                               bigwig.renderer="character",
                               bigwig.sortable="character",
                               bigwig.visible="character",
                               bigwig.windowFunction="character",
                               bigwig.baseline="numeric",
                               bigwig.drawBaseline="character",
                               bigwig.flipAxis="character",
                               bigwig.maximum="numeric",
                               bigwig.minimum="numeric",
                               bigwig.type="character",
                               
                               interval.altColor="character", 
                               interval.color="character",
                               interval.autoScale="character",
                               interval.displayMode="character",                            
                               interval.featureVisibilityWindow="numeric",
                               interval.fontSize="numeric",
                               interval.height="numeric",
                               interval.normalize="character",
                               interval.renderer="character",
                               interval.sortable="character",
                               interval.visible="character",
                               interval.windowFunction="character",
                               
                               bam.altColor="character", 
                               bam.color="character",
                               bam.autoScale="character",
                               bam.displayMode="character",                            
                               bam.featureVisibilityWindow="numeric",
                               bam.fontSize="numeric",
                               bam.showSpliceJunctions="character",
                               bam.colorByTag="character",
                               bam.colorOption="character",
                               bam.flagUnmappedPairs="character",
                               bam.groupByTag="character",
                               bam.maxInsertSize="numeric",
                               bam.minInsertSize="numeric",
                               bam.shadeBasesOption="character",
                               bam.shadeCenters="character",
                               bam.showAllBases="character",
                               bam.sortByTag="character"                            
                     )
)

#' @rdname igvParam
#' @export
igvParam <- function(

  bigwig.altColor="darkgrey", 
  bigwig.color="darkgrey",
  bigwig.autoScale="true",
  bigwig.displayMode="COLLAPSED",
  bigwig.featureVisibilityWindow=-1,
  bigwig.fontSize=10,
  bigwig.normalize="false",
  bigwig.renderer="BAR_CHART",
  bigwig.sortable="true",
  bigwig.visible="true",
  bigwig.windowFunction="mean",
  bigwig.baseline=0.0,
  bigwig.drawBaseline="true",
  bigwig.flipAxis="false",
  bigwig.maximum=50,
  bigwig.minimum=0,
  bigwig.type="LINEAR",
  
  interval.altColor="darkgrey", 
  interval.color="darkgrey",
  interval.autoScale="true",
  interval.displayMode="character",                            
  interval.featureVisibilityWindow=-1,
  interval.fontSize=10,
  interval.height=40,
  interval.normalize="false",
  interval.renderer="BASIC_FEATURE",
  interval.sortable="true",
  interval.visible="true",
  interval.windowFunction="count",
  
  bam.altColor="darkgrey", 
  bam.color="darkgrey",
  bam.autoScale="true",
  bam.displayMode="EXPANDED",                            
  bam.featureVisibilityWindow=-1,
  bam.fontSize=10,
  bam.showSpliceJunctions="false",
  bam.colorByTag="",
  bam.colorOption="UNEXPECTED_PAIR",
  bam.flagUnmappedPairs="false",
  bam.groupByTag="",
  bam.maxInsertSize=1000,
  bam.minInsertSize=50,
  bam.shadeBasesOption="QUALITY",
  bam.shadeCenters="true",
  bam.showAllBases="false",
  bam.sortByTag=""){
  
  igvParamReturn <- new("igvParam",
                        bigwig.altColor=bigwig.altColor,
                        bigwig.color=bigwig.color,
                        bigwig.autoScale=bigwig.autoScale,
                        bigwig.displayMode=bigwig.displayMode,
                        bigwig.featureVisibilityWindow=bigwig.featureVisibilityWindow,
                        bigwig.fontSize=bigwig.fontSize,
                        bigwig.normalize=bigwig.normalize,
                        bigwig.renderer=bigwig.renderer,
                        bigwig.sortable=bigwig.sortable,
                        bigwig.visible=bigwig.visible,
                        bigwig.windowFunction=bigwig.windowFunction,
                        bigwig.baseline=bigwig.baseline,
                        bigwig.drawBaseline=bigwig.drawBaseline,
                        bigwig.flipAxis=bigwig.flipAxis,
                        bigwig.maximum=bigwig.maximum,
                        bigwig.minimum=bigwig.minimum,
                        bigwig.type=bigwig.type,              
                        interval.altColor=interval.altColor, 
                        interval.color=interval.color,
                        interval.autoScale=interval.autoScale,
                        interval.displayMode=interval.displayMode,                            
                        interval.featureVisibilityWindow=interval.featureVisibilityWindow,
                        interval.fontSize=interval.fontSize,
                        interval.height=interval.height,
                        interval.normalize=interval.normalize,
                        interval.renderer=interval.renderer,
                        interval.sortable=interval.sortable,
                        interval.visible=interval.visible,
                        interval.windowFunction=interval.windowFunction,              
                        bam.altColor=bam.altColor, 
                        bam.color=bam.color,
                        bam.autoScale=bam.autoScale,
                        bam.displayMode=bam.displayMode,                            
                        bam.featureVisibilityWindow=bam.featureVisibilityWindow,
                        bam.fontSize=bam.fontSize,
                        bam.showSpliceJunctions=bam.showSpliceJunctions,
                        bam.colorByTag=bam.colorByTag,
                        bam.colorOption=bam.colorOption,
                        bam.flagUnmappedPairs=bam.flagUnmappedPairs,
                        bam.groupByTag=bam.groupByTag,
                        bam.maxInsertSize=bam.maxInsertSize,
                        bam.minInsertSize=bam.minInsertSize,
                        bam.shadeBasesOption=bam.shadeBasesOption,
                        bam.shadeCenters=bam.shadeCenters,
                        bam.sortByTag=bam.sortByTag)                                         
  
  
  return(igvParamReturn)
}


implode <- function (strings, sep = "", finalSep = NULL, ranges = FALSE)
{
  # Transform runs of integers into ranges
  # This is surprisingly tricky to get right!
  if (ranges && is.integer(strings))
  {
    # Perform run-length encoding on the differences between elements
    gapRunLengths <- rle(diff(strings))
    
    # Mark all elements not taken and find ranges (>1 consecutive unit difference)
    taken <- rep(FALSE, length(strings))
    withinRange <- gapRunLengths$values == 1 & gapRunLengths$lengths > 1
    
    # Convert range groups into strings, marking elements as taken to avoid double-counting
    rangeStrings <- lapply(which(withinRange), function(i) {
      # NB: Sum of a length-zero vector is zero
      start <- sum(gapRunLengths$lengths[seq_len(i-1)]) + 1
      end <- start + gapRunLengths$lengths[i]
      taken[start:end] <<- TRUE
      return (paste(strings[start], strings[end], sep="-"))
    })
    
    # Convert remaining elements into strings
    nonRangeStrings <- lapply(which(!withinRange), function(i) {
      start <- sum(gapRunLengths$lengths[seq_len(i-1)]) + 1
      end <- start + gapRunLengths$lengths[i]
      toKeep <- setdiff(start:end, which(taken))
      taken[toKeep] <<- TRUE
      return (as.character(strings)[toKeep])
    })
    
    # Arrange list of strings in the right order, and convert back to character vector
    strings <- vector("list", length(withinRange))
    strings[withinRange] <- rangeStrings
    strings[!withinRange] <- nonRangeStrings
    strings <- unlist(strings)
  }
  else
    strings <- as.character(strings)
  
  if (length(strings) == 1)
    return (strings[1])
  else if (length(strings) > 1)
  {
    result <- strings[1]
    for (i in 2:length(strings))
    {
      if (i == length(strings) && !is.null(finalSep))
        result <- paste(result, strings[i], sep=finalSep)
      else
        result <- paste(result, strings[i], sep=sep)
    }
    return (result)
  }
}


relativePath <- function (path, referencePath)
{
  mainPieces <- strsplit(expandFileName(path), .Platform$file.sep, fixed=TRUE)[[1]]
  refPieces <- strsplit(expandFileName(referencePath), .Platform$file.sep, fixed=TRUE)[[1]]
  
  shorterLength <- min(length(mainPieces), length(refPieces))
  firstDifferentPiece <- min(which(mainPieces[1:shorterLength] != refPieces[1:shorterLength])[1], shorterLength, na.rm=TRUE)
  newPieces <- c(rep("..", length(refPieces)-firstDifferentPiece), mainPieces[firstDifferentPiece:length(mainPieces)])
  
  return (implode(newPieces, sep=.Platform$file.sep))
}

expandFileName <- function (fileName, base = getwd())
{
  fileName <- path.expand(fileName)
  
  # A leading slash, with (Windows) or without (Unix) a letter and colon, indicates an absolute path
  fileName <- ifelse(fileName %~% "^([A-Za-z]:)?/", fileName, file.path(base,fileName))
  
  # Remove all instances of '/.' (which are redundant), recursively collapse
  # instances of '/..', and remove trailing slashes
  fileName <- gsub("/\\.(?=/)", "", fileName, perl=TRUE)
  while (length(grep("/../", fileName, fixed=TRUE) > 0))
    fileName <- sub("/[^/]*[^./][^/]*/\\.\\.(?=/)", "", fileName, perl=TRUE)
  if (length(grep("/..$", fileName, perl=TRUE) > 0))
    fileName <- sub("/[^/]*[^./][^/]*/\\.\\.$", "", fileName, perl=TRUE)
  fileName <- gsub("/*\\.?$", "", fileName, perl=TRUE)
  
  return (fileName)
}


