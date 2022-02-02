
def test_coordination_T(example_gene_coor,res_coor):
    res=example_gene_coor.gene_coordination_test(test="fisher")
    assert all([res[i]==res_coor[i] for i in range(2,13)])

def test_coordination_F(example_gene_uncoor,res_uncoor):
    res=example_gene_uncoor.gene_coordination_test(test="fisher")
    assert all([res[i]==res_uncoor[i] for i in range(2,13)])