import os.path

from matplotlib import pyplot as plt
from matplotlib import image as img


class PlotMultiple():

    def __init__(self):
        self.ax = None
        self.fig = None
        self.index_row = 0
        self.index_col = 0
        self.nrow = 1
        self.ncol = 1
        self.axhere = None

    def start_multiple_plot(self, nrow, ncol):
        self.nrow = nrow
        self.ncol = ncol
        # self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(8, 6), frameon=False,
        #                                  gridspec_kw={'wspace': 0, 'hspace': 0})

        # self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(15.9, 18), frameon=False,gridspec_kw={'wspace': 0, 'hspace': 0})

        # self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(16, 12), frameon=False,gridspec_kw={'wspace': 0, 'hspace': 0})
        self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(16, 12), frameon=False,
                                         gridspec_kw={'wspace': 0, 'hspace': 0})

        # self.ax.set_axis_off()
        # self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(8, 6), frameon=False,
        #                                  gridspec_kw={'wspace': 0, 'hspace': 0})

    def start_multiple_plot_advanced(self, nrow, ncol, xfigsize, yfigsize, wspace, hspace, frameon):
        self.nrow = nrow
        self.ncol = ncol
        self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(xfigsize, yfigsize), frameon=frameon,
                                         gridspec_kw={'wspace': wspace, 'hspace': hspace})

    def start_multiple_plot_polar(self, nrow, ncol, xfigsize, yfigsize, wspace, hspace, frameon):
        self.nrow = nrow
        self.ncol = ncol
        self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(xfigsize, yfigsize), frameon=frameon,
                                         gridspec_kw={'wspace': wspace, 'hspace': hspace},
                                         subplot_kw={'projection': 'polar'})

    def set_text(self, x, y, s):
        plt.text(x, y, s, fontsize=10, backgroundcolor='w')


    def set_text_size(self, x, y, s,fontsize):
        plt.text(x, y, s, fontsize=fontsize, backgroundcolor='w')

    def get_axes(self, index_row, index_col):
        if self.nrow > 1 and self.ncol > 1:
            axhere = self.ax[index_row, index_col]
        elif self.nrow == 1 and self.ncol > 1:
            axhere = self.ax[index_col]
        if self.nrow > 1 and self.ncol == 1:
            axhere = self.ax[index_row]
        return axhere

    def plot_blank(self, index_row, index_col):
        if self.nrow > 1 and self.ncol > 1:
            axhere = self.ax[index_row, index_col]
        elif self.nrow == 1 and self.ncol > 1:
            axhere = self.ax[index_col]
        if self.nrow > 1 and self.ncol == 1:
            axhere = self.ax[index_row]
        axhere.axis('off')

    def plot_blank_with_title(self, index_row, index_col, title):
        if self.nrow > 1 and self.ncol > 1:
            axhere = self.ax[index_row, index_col]
        elif self.nrow == 1 and self.ncol > 1:
            axhere = self.ax[index_col]
        if self.nrow > 1 and self.ncol == 1:
            axhere = self.ax[index_row]

        axhere.set_xticks([])
        axhere.set_yticks([])

        if title is not None:
            axhere.set_title(title)

    def plot_image_hypernets(self,file_img,index_row,index_col,title):

        if self.nrow > 1 and self.ncol > 1:
            axhere = self.ax[index_row, index_col]
        elif self.nrow == 1 and self.ncol > 1:
            axhere = self.ax[index_col]
        if self.nrow > 1 and self.ncol == 1:
            axhere = self.ax[index_row]


        from PIL import Image
        image = Image.open(file_img)
        rimage = image.rotate(270,expand=True)

        axhere.imshow(rimage)

        w, h = rimage.size

        ##central point
        axhere.axvline(w / 2, 0.48,0.52,color='red', linewidth=0.25)
        axhere.axhline(h / 2, 0.48, 0.52, color='red', linewidth=0.25)
        ##grid
        incremx = int(w/4)
        incremy = int(h/4)
        for x in range(0,w,incremx):
            axhere.axvline(x, color='red', linewidth=0.5)
        for y in range(0,h,incremy):
            axhere.axhline(y, color='red', linewidth=0.5)

        if title is not None:
            axhere.set_title(title)
        axhere.set_xticks([])
        axhere.set_yticks([])


    def plot_image_title(self, file_img, index_row, index_col, title):
        image = img.imread(file_img)
        if self.nrow > 1 and self.ncol > 1:
            axhere = self.ax[index_row, index_col]
        elif self.nrow == 1 and self.ncol > 1:
            axhere = self.ax[index_col]
        if self.nrow > 1 and self.ncol == 1:
            axhere = self.ax[index_row]

        image_end = image

        axhere.imshow(image_end)
        if title is not None:
            axhere.set_title(title)
        axhere.set_xticks([])
        axhere.set_yticks([])

    def plot_image(self, file_img, index_row, index_col):
        image = img.imread(file_img)
        if self.nrow > 1 and self.ncol > 1:
            axhere = self.ax[index_row, index_col]
        elif self.nrow == 1 and self.ncol > 1:
            axhere = self.ax[index_col]
        if self.nrow > 1 and self.ncol == 1:
            axhere = self.ax[index_row]

        image_end = image

        axhere.imshow(image_end)

        axhere.axis(False)

    def tigth_layout(self):
        self.fig.tight_layout()

    def set_global_legend(self, handles, str_legend):
        # self.fig.legend(handles, str_legend, fontsize = 8, loc='lower center', ncol=len(str_legend), markerscale=1.5,bbox_to_anchor=(0.55,0.08))

        self.fig.legend(handles, str_legend, fontsize=8, loc='lower center', ncol=len(str_legend), markerscale=1.0,
                        bbox_to_anchor=(0.55, 0.06))

        # self.fig.legend(handles, str_legend, fontsize=8, loc='lower center', ncol=len(str_legend), markerscale=1.0)

    def set_global_legend_2(self, handles, str_legend):
        # self.fig.legend(handles, str_legend, fontsize = 8, loc='lower center', ncol=len(str_legend), markerscale=1.5,bbox_to_anchor=(0.55,0.08))
        self.fig.legend(handles, str_legend, fontsize=11, loc='lower center', ncol=len(str_legend), markerscale=1.0,
                        bbox_to_anchor=(0.50, -0.03))

    def set_global_legend_3(self, handles, str_legend):
        # self.fig.legend(handles, str_legend, fontsize = 8, loc='lower center', ncol=len(str_legend), markerscale=1.5,bbox_to_anchor=(0.55,0.08))
        handles_new = [handles[0], handles[4], handles[1], handles[5], handles[2], handles[6], handles[3], handles[7]]
        str_legend_new = [str_legend[0], str_legend[4], str_legend[1], str_legend[5], str_legend[2], str_legend[6],
                          str_legend[3], str_legend[7]]

        self.fig.legend(handles_new, str_legend_new, fontsize=11, loc='lower center', ncol=4, markerscale=1.0,
                        bbox_to_anchor=(0.53, 0.04), frameon=False)

    def save_fig(self, file_out):
        if file_out.endswith('.tif'):
            # print('es aqui')
            plt.savefig(file_out, dpi=300, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})
            # plt.savefig(file_out, dpi=300,bbox_inches='tight',transparency=False,facecolor='white',pil_kwargs={"compression": "tiff_lzw"})
            # plt.savefig(file_out, dpi=300, bbox_inches='tight', transparency=False, facecolor='white',pil_kwargs={"compression": "tiff_lzw"})
        else:
            # plt.savefig(file_out, dpi=300, bbox_inches='tight', transparency=False, facecolor='white')
            plt.savefig(file_out, dpi=300, bbox_inches='tight', facecolor='white')
        # plt.savefig(file_out, dpi=300)



    def close_plot(self):
        plt.close()
