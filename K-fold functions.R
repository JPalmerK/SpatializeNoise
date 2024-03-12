

####################################################
## Create the plot for signal excess##
###################################################
figFilesloc ='C:\\BitBucketRepositories\\CABOWProcessing\\R\\Publication Figures'

My_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 12),
  legend.title=element_text(size=10))
# 
# # Example of how to set up figures
# ragg::agg_jpeg(paste(figFilesloc,"\\ Figure 3 NL.jpg",sep=""), 
#                width = 25, 
#                height = 16, 
#                units = "cm", res = 300)
# NLPlot
# dev.off()


######################################################################
# Build and cross validate models
#####################################################################

fieldModelResults <-function(NLData = simplifNL,
                             band = '500Hz', 
                             predDrift ='ADRIFT_046'){
  # Function for pulling model and prediction data from kfold
  
   modelData = subset(NLData, band = '500Hz', DriftName !=predDrift)
   predData = subset(NLData, band = '500Hz', DriftName ==predDrift)
  
   model.loc<-cbind(modelData$Lon,
                  modelData$Lat) #locations
   pred.loc<-cbind(predData$Lon,
                   predData$Lat) #locations
   
   model.y = modelData$RelNoise[modelData$DriftName !=predDrift]
   pred.y = predData$RelNoise[predData$DriftName !=predDrift]
   
  
   out<-list(model.loc, model.y, pred.loc, pred.y)
  return(out)
  
}