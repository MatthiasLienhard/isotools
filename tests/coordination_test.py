
def test_coordination_T(example_gene_coor, res_coor):
    res = example_gene_coor.gene_coordination_test(test="fisher")
    assert all([res[0][i] == res_coor[0][i] for i in range(2, 14)])
    assert res[0][0] > 0 and res[0][0] < .05


def test_coordination_F(example_gene_uncoor, res_uncoor):
    res = example_gene_uncoor.gene_coordination_test(test="fisher")
    assert all([res[0][i] == res_uncoor[0][i] for i in range(2, 14)])
    assert res[0][0] > .3 and res[0][0] < .5