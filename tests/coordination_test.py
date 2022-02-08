import numpy as np

def test_coordination(example_gene_coor):
    example_gene_coor.data['coverage'] = np.array([300, 60, 100, 350])
    res_pos = example_gene_coor.gene_coordination_test(test="fisher")
    assert res_pos[0][0] >= 0 and res_pos[0][0] < .05 , 'Test should yield significant p-value'
    example_gene_coor.data['coverage'] = np.array([310, 380, 310, 350])
    res_neg = example_gene_coor.gene_coordination_test(test="fisher")
    assert res_neg[0][0] > .1 and res_neg[0][0] <= 1, 'Test should not yield significant p-value'