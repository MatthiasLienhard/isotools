import numpy as np

def test_coordination(example_gene_coor):
    example_gene_coor.data['coverage'] = np.array([[300, 60, 100, 350]])
    res_pos = example_gene_coor.coordination_test(test="fisher")
    assert res_pos[0][8] >= 0 and res_pos[0][8] < .05 , 'Test should yield significant p-value'
    example_gene_coor.data['coverage'] = np.array([[310, 380, 310, 350]])
    res_neg = example_gene_coor.coordination_test(test="fisher")
    assert res_neg[0][8] > .1 and res_neg[0][8] <= 1, 'Test should not yield significant p-value'