from model import *

# Data visualization
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

class Show(Model):
    def __init__(self):
        super().__init__()

    def show_positions(self, seq, top, dic):
        '''
        This function allows to plot the results. The plot is saved in the /output/prediction folder.

        :param seq: wether the resilts are from a
        :param top: The top n predictions to display
        :param dic: a dic in format {gene name: [start position, end position, strand (p or q),[intron position]}
        :return: output the plot the plot and
        '''
        # dic = {gene name: [start,end,strand]}
        # If the seq argument is setted to fasta
        if seq == "fasta":

            # We create a figure
            f = plt.figure(figsize=(30, 22))

            # And an axis
            ax = f.add_subplot(311)

            # Then we sort predictions lists according to the predicted probability, we also add their index positions
            self.p_preds = sorted([(self.p_preds[i], i) for i in range(len(self.p_preds))], key=lambda x: x[0])
            self.q_preds = sorted([(self.q_preds[i], i) for i in range(len(self.q_preds))], key=lambda x: x[0])

            # Finally, we will only draw the top genes (top passed as arguments) for both strands
            p_top = self.p_preds[len(self.p_preds) - top:len(self.q_preds)]
            q_top = self.q_preds[len(self.q_preds) - top:len(self.q_preds)]
            # Then we plot the data
            ax.plot([2 for i in self.p_preds], linewidth=15, color="blue", label="Direct 5' -> 3'")
            ax.plot([-2 for i in self.q_preds], linewidth=15, color="grey", label="Reverse 5' -> 3'")
            ax.set_yticks((3,2,-2,-3))
            ax.set_yticklabels(("p strand","5'->3'","q_strand"," 5'->3'"))
            ax.set_ylim((-6,6))
            ax.legend()
            # For each probabilities
            for i in p_top:
                # We draw a line
                ax.vlines(ymin=1.8, ymax=2.2, x=i[1], colors='red', linewidth=2)
            # Same for the q strand
            for i in q_top:
                ax.vlines(ymin=-1.8, ymax=-2.2, x=i[1], colors='red', linewidth=2)

            # if dic is passed as an argument, we plot the gene positions
            if isinstance(dic, dict):
                # We plot the predicted sites
                for i in dic:
                    if dic[i][2] == 'p':
                        ymin = 1
                        ymax = 3
                        start = dic[i][0]
                        end = dic[i][1]
                        st = 2
                    elif dic[i][2] == 'q':
                        ymin = -1
                        ymax = -3
                        start = dic[i][1]
                        end = dic[i][0]
                        st = -2
                    else:
                        print(
                            "Format error, please make sure the gene annotations are provided in a dict format such as = {gene name:[start,end,strand]}")
                    # We plot the TSS of the provided genes
                    ax.vlines(ymin=ymin, ymax=ymax, x=start, colors='green', linewidth=6)
                    ax.text(start + 15, ymax -0.5 , "{}".format(i), fontsize='x-large')
                    # Then its end
                    ax.vlines(ymin=ymin, ymax=ymax, x=end, colors='orange', linewidth=6)
                    ax.text(end + 15, ymax - 0.5, "{}".format(i), fontsize='x-large')
                    for j in dic[i][3]:
                        ax.plot([m for m in j],[st for u in range(len(j))], linewidth=15, color="red")

            # We set the axis limit such as the results are visibles
            ax.set_ybound(-4, 4)
            # Finally, we set the legend
            legend_elements = [Line2D([0], [0], color='r', lw=2, label='Predicted {} sites'.format(self.motif_names[0])),
                               Line2D([0], [0], color='g', label='TSS'), Line2D([0], [0], color='orange', label='STOP'),
                               Patch(facecolor='red', edgecolor='r',
                                     label='Intron')
                               ]
            ax.legend(handles=legend_elements, loc='upper right')

            ax1 = f.add_subplot(312)
            ax2 = f.add_subplot(313)

            clist = [(0.0, "blue"), (1, "red")]

            self.p_preds = sorted(self.p_preds, key=lambda x: x[0][0])
            self.p_preds = self.p_preds[len(self.p_preds) - top:len(self.p_preds)]
            self.q_preds = sorted(self.q_preds, key=lambda x: x[0][0])
            self.q_preds = self.q_preds[len(self.q_preds) - top:len(self.q_preds)]

            x1, y1 = [i[1] for i in self.p_preds], np.array([i[0][0] for i in self.p_preds])
            x2, y2 = [i[1] for i in self.q_preds], np.array([i[0][0] for i in self.q_preds])
            print("Strict motif:")
            print(DataLoader.get_strict_motif(self))
            DataLoader.normalize_matrix(self)
            print("Predicted motif for the direct strand")
            for i in range(len(x1)):
                s = self.p_seq[x1[i]:x1[i] + len(self.A)]
                print("S = {}, score = {}".format(super().reverse(s),
                                                  super().get_alignment_score(s)))
                y1[i]=y1[i]*DataLoader.get_alignment_score(self, s)
            for i in range(len(x2)):
                s = self.p_seq[x2[i]:x2[i] + len(self.A)]
                print("S = {}, score = {}".format(super().reverse(s),
                                                  super().get_alignment_score(s)))
                y2[i]=y2[i]*DataLoader.get_alignment_score(self, s)


            from matplotlib import cm
            ax1 = ax1.scatter(x1, y1, c=y1, cmap='RdBu')
            ax2 = ax2.scatter(x2, y2, c=y2, cmap='RdBu')

            f.colorbar(ax1,orientation="horizontal")
            plt.savefig(os.path.join(self.prediction_dir, " {}.png".format([i for i in dic][0])))
            plt.show()


            # ax1.bar(x1,y1,cmap=plt.get_cmap('YlOrBr'))
            # ax1.scatter([i[1] for i in self.p_preds],[i[0] for i in self.p_preds])

            # ax2.bar(x2,y2,color= rvb(y2))
            # ax2.scatter([i[1] for i in self.q_preds],[np.log(np.arcsin(i[0])+1) for i in self.q_preds])

class Experiment(Show):
    def __init__(self):
        super().__init__()
