import pickle
import soa_motif_analysis


#test_motif = soa_motif_analysis.get_motifs(meme_data_dir='/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/issac-workspace-IE/soa/meme_bin/CLSTR1667_meme_out/', e_val_threshold=1)

#pickle.dump(test_motif, open('test.p', 'wb'))

test_motif = pickle.load(open('test.p', 'rb'))

print(soa_motif_analysis.calculate_motif_distance(test_motif[0], test_motif[0], padded=False))
