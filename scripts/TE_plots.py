#!/bin/python
# This script has modules for plotting TE anootations mapped to a consensus sequence. 
#
#
#
# Abin Abraham
# created on: 2017-12-20 19:18:182


def plot_MapToConsensus(consArray, element_NAME, annotation_NAME, outputDir, saveFlag):

    import matplotlib.pyplot as plt
    #from cycler import cycler
    consensusLength= consArray.shape[1]
    consBases = np.arange(0.0, consensusLength, 1, dtype=int)
    meanConsArray = np.mean(consArray, axis=0)
    stdConsArray = np.std(consArray, axis=0)

    plt.figure(figsize=(8, 5), dpi=160)
    
    plt.subplot(2, 1, 1) #RAW DATA l
    #plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))
    for i in range(len(consArray)):
        plt.plot(consBases, consArray[i], 'b-', linewidth = 1, alpha=0.1)
    plt.grid()
    plt.title('Genomic Instances of ' + element_NAME + ' Annotated \n with ' + annotation_NAME )
    plt.ylabel(annotation_NAME)
    [x1,x2,y1,y2] = plt.axis()
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])

    plt.subplot(2, 1, 2) #MEAN and STD 
    plt.plot(consBases, meanConsArray, 'r-', linewidth = 1.2, alpha=1)
    plt.plot(consBases, meanConsArray + 2*stdConsArray, 'b:', linewidth = 0.7, alpha=0.8)
    plt.plot(consBases, meanConsArray - 2*stdConsArray, 'b:', linewidth = 0.7, alpha=0.8)
    plt.axis((x1,x2,-1*y2,y2))
    plt.grid()
    plt.title('Mean (Red) and 2*SD (blue) of ' + annotation_NAME + ' annotation of ' + element_NAME)
    plt.xlabel('consensus bases (1start)')
    plt.ylabel(annotation_NAME)
    saveFigName = os.path.join(outputDir, element_NAME+"_"+annotation_NAME+"_MappedToConsensus_" + str(datetime.date.today())+ ".png")
    plt.subplots_adjust(hspace=0.5)
    plt.show()
    if saveFlag: 
        plt.savefig(saveFigName)

    plt.close()