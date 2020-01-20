# ------------------------------------------------
#                 GRAPH MODULE
# ------------------------------------------------
# Function to plot Blastp results.

# MODULE IMPORTATION
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

# FUNCTION DEFINITION
def blastp(output,nseqs,runid):
    """output: name of the tsv file with the Blastp results (str).
    nseqs: number of query sequences (int).
    runid: numeric id of the program's run (int).
    
    This funcion creates a graph for each query sequence and its
    Blastp hits, showing query coverage, identity % and evalue.
    """
    
    # Repeat process for each query.
    for i in range(0,nseqs):
        qseq = []
        query = False
        j = 0
        with open(output,"r") as results:
            for line in results:
                data = line.split("\t")
                if line[0]!="#":
                    if data[0] not in qseq:
                        if j == i:
                            qseq.append(data[0])
                            query = True
                            sseq,qcov,pident,evalue = [],[],[],[]
                        else:
                            qseq.append(data[0])
                            query = False
                        j += 1
                    
                    if query == True:
                        sseq.append(data[1])
                        qcov.append(data[2])
                        pident.append(data[3])
                        evalue.append(data[4])

        qcov = np.asarray(qcov,dtype=int)
        plot_filename = "./results/blastp/blastp" + str(runid) + ".graph" + str(i) + ".png"
        y_pos = np.arange(len(sseq))

        # RGB color gradient to represent evalue.
        gradient = [(171,235,118),(130,224,170),(88,214,141),(46,204,113),(40,180,99),(35,155,86)]
        gradient = np.divide(gradient,255)
        evalue = np.asarray(evalue,dtype=float)

        # Assign bar color according to e-value of the hit
        color = []
        for x in range(0,len(sseq)):
            if evalue[x] < 1e-80:
                color.append(gradient[5])
            elif (evalue[x] > 1e-80) and (evalue[x] < 1e-60):
                color.append(gradient[4])
            elif (evalue[x] > 1e-60) and (evalue[x] < 1e-40):
                color.append(gradient[3])
            elif (evalue[x] > 1e-40) and (evalue[x] < 1e-20):
                color.append(gradient[2])
            elif (evalue[x] > 1e-20) and (evalue[x] < 1):
                color.append(gradient[1])
            else:
                color.append(gradient[0])

        # Make bar plot of query coverage values.
        plt.barh(y_pos,qcov,align='center',color=color)

        # Remove top and right axes.
        ax=plt.gca()
        for pos in ['right','top']:
            ax.spines[pos].set_visible(False)

        # e-value color legend.
        ev = mpatches.Patch(color="white", label='e-value')
        green0 = mpatches.Patch(color=gradient[0], label='> 1')
        green1 = mpatches.Patch(color=gradient[1], label='(1e-20;1)')
        green2 = mpatches.Patch(color=gradient[2], label='(1e-40;1e-20)')
        green3 = mpatches.Patch(color=gradient[3], label='(1e-60;1e-40)')
        green4 = mpatches.Patch(color=gradient[4], label='(1e-80;1e-60)')
        green5 = mpatches.Patch(color=gradient[5], label='< 1e-80')
        plt.legend(handles=[ev,green0,green1,green2,green3,green4,green5],bbox_to_anchor=(1,1), loc='upper left', borderaxespad=0.,fontsize=8)

        # Show identity % labels.
        for j in range(0,len(pident)):
            pos = np.array((qcov[j],j))
            name = "  " + str(pident[j]) + "% id"
            plt.text(pos[0],pos[1],name,fontsize=7)

        # Titles and y axis names.
        plt.yticks(y_pos,sseq)
        plt.xlabel('Query coverage')
        plt.ylabel('Subject sequence')
        plt.title('BLASTP HITS. Query seq: %s\n' % qseq[i],fontweight="bold")

        # Save png image of each plot.
        plt.savefig(plot_filename,bbox_inches='tight')
        # Close figure to make a new one.
        plt.clf()
