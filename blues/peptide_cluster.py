import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy
from scipy.spatial.distance import squareform


def main():
    # input files
    traj_file = 'coupled.xtc'
    gro_file = 'coupled.gro'

    # load the entire trajectory
    trj = None
    for chunk in md.iterload(traj_file, top=gro_file, chunk=500, stride=1):
        temp_trjs = chunk
        temp_trjs.remove_solvent(inplace=True)
        if trj is None:
            trj = temp_trjs
        else:
            trj = trj + temp_trjs
        del temp_trjs

    # atoms for clustering
    index_rmsd = trj.topology.select('name CA')

    # superpose
    trj.superpose(reference=trj[0], atom_indices=index_rmsd, ref_atom_indices=index_rmsd)

    # calculating the pair-wise rmsd of each frame to be used for clustering
    distances = np.empty((trj.n_frames, trj.n_frames))
    for i in range(trj.n_frames):
        distances[i] = md.rmsd(trj, trj, i, atom_indices=index_rmsd)

    reduced_distances = squareform(distances, checks=False)
    linkage = scipy.cluster.hierarchy.linkage(reduced_distances, method='average')

    plt.title('Interaction Distance Average linkage hierarchical clustering')
    scipy.cluster.hierarchy.dendrogram(linkage, no_labels=True, count_sort='descendent')
    plt.show()
    k = int(input('Enter k value (max number of clusters): '))
    crit = 'maxclust'
    clusterindices = scipy.cluster.hierarchy.fcluster(linkage, k, criterion=crit)

    # save pdb files for each cluster
    index_count = range(0, len(clusterindices))
    for ii in range(1, max(clusterindices) + 1):
        cluster = [i for i in index_count if clusterindices[i] == ii]
        trj[cluster].save_pdb('cluster_' + str(ii) + '.pdb')


if __name__ == "__main__":
    main()
