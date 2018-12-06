from __future__ import division
import re
import ete3

class phydelity_output(object):
    """
    Print output cluster and tree files
    """

    def __init__(self, ori_tree_string=None, leaf_node_id_to_leafname=None, taxon_to_clusterid=None, clusterid_to_taxa=None, outfname=None):

        self.ori_tree_string = ori_tree_string

        # list of leaves by name
        self.taxon_list = leaf_node_id_to_leafname.values()

        # convert leaf node_id to leaf name
        self.taxon_to_clusterid = {leaf_node_id_to_leafname[taxon]:clusterid for taxon, clusterid in taxon_to_clusterid.items()}
        self.clusterid_to_taxa = {clusterid:[leaf_node_id_to_leafname[taxon] for taxon in taxa] for clusterid, taxa in clusterid_to_taxa.items()}

        self.outfname = outfname

    def cluster_output(self):

        curr_tree_string = self.ori_tree_string

        # get node-ids annotated for internal nodes first
        edited_tree_string = []
        prev_end = 0
        for expr in re.finditer('\)(\d+):', curr_tree_string):
            edited_tree_string.append(curr_tree_string[prev_end:expr.start()+1])
            edited_tree_string.append('[&NODE_ID=\'{}\']:'.format(expr.group(1)))
            prev_end = expr.end()
        edited_tree_string.append(curr_tree_string[prev_end:])
        curr_tree_string = ''.join(edited_tree_string)

        with open('cluster_{}.txt'.format(self.outfname), 'w') as output:

            output.write('CLUSTER\tTAXA\r\n')
            for taxon, clusterid in self.taxon_to_clusterid.items():

                output.write('{}\t{}\r\n'.format(clusterid, re.sub("(^'|'$)", '', taxon)))

                taxon_start_index = curr_tree_string.find(taxon)
                if taxon_start_index < 0:
                    raise SystemExit('\r\nERROR: Problem parsing taxon ({}) in tree string.\r\n'.format(taxon))
                else:
                    curr_tree_string = curr_tree_string.replace(taxon, "'{}'[&CLUSTER_ID='{}',CLUSTERED=1]".format(re.sub("(^'|'$)", "", taxon), clusterid))

            for taxon in self.taxon_list:
                if taxon not in self.taxon_to_clusterid:
                    curr_tree_string = curr_tree_string.replace(taxon, "'{}'[&CLUSTERED=0]".format(re.sub("(^'|'$)", "", taxon), clusterid))

        return curr_tree_string

    def figtree_output(self, modified_tree_string):

        with open('tree_{}.tre'.format(self.outfname), 'w') as output:
            output.write('#NEXUS\r\nBegin taxon;\r\n\tDimensions ntax={};\r\n\t\tTaxlabels\r\n'.format(len(self.taxon_list)))
            for taxon in self.taxon_list:
                if taxon in self.taxon_to_clusterid:
                    output.write("\t\t\t'{}_cluster{}'\r\n".format(re.sub("(^'|'$)", "", taxon), self.taxon_to_clusterid[taxon]))
                else:
                    output.write('\t\t\t{}\r\n'.format(taxon))

            output.write('\t\t\t;\r\nEnd;\r\nBegin trees;\r\n\tTranslate\r\n')

            for i, taxon in enumerate(self.taxon_list):
                if taxon in self.taxon_to_clusterid:
                    output.write("\t\t{:>4} '{}_cluster{}'{}\r\n".format(i+1, re.sub("(^'|'$)", '', taxon), self.taxon_to_clusterid[taxon], '' if i+1 == len(self.taxon_list) else ','))
                else:
                    output.write("\t\t{:>4} '{}'{}\r\n".format(i+1, re.sub("(^'|'$)", '', taxon), '' if i+1 == len(self.taxon_list) else ','))

                taxon_start_index = modified_tree_string.find("'{}'".format(re.sub("(^'|'$)", '', taxon)))
                if taxon_start_index < 0:
                    raise SystemExit('\r\nERROR: Problem parsing taxon ({}) in tree string.\r\n'.format(taxon))
                else:
                    modified_tree_string = modified_tree_string.replace("'{}'".format(re.sub("(^'|'$)", "", taxon)), str(i+1))

            output.write(';\r\ntree TREE1 = {}\r\nEnd;\r\n'.format(modified_tree_string))

    def generate_heatmap(self, ete3_tree, *args):

        for node in ete3_tree:
            if node.is_leaf():
                taxon = node.name

                tc_index = 0
                for args_index in range(0, len(args), 2):
                    taxon_to_classification = args[args_index]
                    class_to_color = args[args_index+1]

                    try:
                        classification = taxon_to_classification[taxon]
                        color = class_to_color[classification]
                    except:
                        classification = ''
                        color = 'white'

                    class_face = ete3.TextFace(classification, ftype='Arial', fsize=2, fgcolor='white')
                    class_face.background.color = color
                    class_face.hz_align = 1
                    class_face.vt_align = 1
                    class_face.margin_left = 5
                    class_face.margin_right = 5
                    node.add_face(class_face, tc_index, "aligned")
                    tc_index += 1

        return ete3_tree

    def generate_color_scheme(self, class_list):
        from colorsys import hls_to_rgb
        import random
        random.seed(666) # maintain consistent h_list shuffle

        def hls2hex(h, l, s):
            ''' Converts from HLS to RGB '''
            return '#%02x%02x%02x' %tuple(map(lambda x: int(x*255), hls_to_rgb(h, l, s)))
        def get_color(h, s=0.5, l=0.5):
            ''' Returns a RGB color code with the specific Hue, Saturation and
            lightness.'''
            return hls2hex(h, l, s)

        # get widest increment of color hues
        h_increment = 1/(len(class_list))
        h_list = [0+(cl*h_increment) for cl in range(len(class_list))]
        random.shuffle(h_list)
        return {classification:get_color(h_list[cl]) for cl, classification in enumerate(class_list)}

    def ete3_pdf_tree_output(self):
        output_tree = ete3.Tree(self.ori_tree_string)
        output_tree.ladderize()

        # treestyle
        ts = ete3.TreeStyle()
        ts.show_leaf_name = False

        # generate color scheme for taxon_to_clusterid
        if len(self.clusterid_to_taxa.keys()) > 1:
            clusterid_to_color = self.generate_color_scheme(self.clusterid_to_taxa.keys())
        else:
            clusterid_to_color = {self.clusterid_to_taxa.keys()[0]:'#ff2929'}

        for n, node in enumerate(output_tree.traverse(strategy='levelorder')):
            if n == 0:
                try:
                    ts.scale_length = float('{:.3f}'.format(node.get_farthest_leaf()[-1]/10))
                except:
                    pass

            if node.is_leaf():
                # color branches
                ns = ete3.NodeStyle()
                ns["size"] = 0 # no node shape

                taxon = node.name
                if taxon in self.taxon_to_clusterid:
                    clusterid = self.taxon_to_clusterid[taxon]
                    ns["hz_line_color"] = clusterid_to_color[clusterid]
                    # write taxon names aligned to the right
                    taxon_name = ete3.TextFace(taxon, ftype='Arial', fsize=2, bold=True, fgcolor=clusterid_to_color[clusterid])
                else:
                    # write taxon names aligned to the right
                    taxon_name = ete3.TextFace(taxon, ftype='Arial', fsize=2, fstyle="italic")

                node.set_style(ns)
                taxon_name.margin_left = 2
                node.add_face(taxon_name, column=0, position='branch-right')

            else:
                ns = ete3.NodeStyle()
                ns["size"] = 0 # no node shape

                # set node style
                node.set_style(ns)

        heatmap_headers = ['Cluster-ID']
        output_tree = self.generate_heatmap(output_tree, self.taxon_to_clusterid, clusterid_to_color)

        # heatmap header
        for lh_index, legend_header in enumerate(heatmap_headers):
            header_face = ete3.TextFace(legend_header, ftype='Arial', fsize=2)
            header_face.hz_align = 1
            header_face.vt_align = 1
            header_face.margin_left = 5
            header_face.margin_right = 5
            ts.aligned_header.add_face(header_face, lh_index)

        # render as pdf
        output_tree.render('pdftree_{}.pdf'.format(self.outfname), tree_style=ts)