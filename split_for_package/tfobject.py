class TFObj:
    def __init__(self,
    tf_activities_condition : list,
    tf_activities_cluster : list,
    average_gene_expression : list,
    regulon : pd.DataFrame,
    CTR_input_condition : list,
    CTR_input_cluster : list,
    intracellular_network_condition : list,
    intracellular_network_cluster : list):

        self.tf_activities_condition = tf_activities_condition
        self.tf_activities_cluster = tf_activities_cluster
        self.average_gene_expression = average_gene_expression
        self.regulon = regulon
        self.CTR_input_condition = CTR_input_condition
        self.CTR_input_cluster = CTR_input_cluster
        self.intracellular_network_condition = intracellular_network_condition
        self.intracellular_network_cluster = intracellular_network_cluster

def make_TFOBj(tf_activities_condition : list,
    tf_activities_cluster : list,
    average_gene_expression : list,
    regulon : pd.DataFrame,
    CTR_input_condition : list,
    CTR_input_cluster : list,
    intracellular_network_condition : list,
    intracellular_network_cluster : list):

        tf = TFObj(tf_activities_condition,
            tf_activities_cluster,
            average_gene_expression,
            regulon,
            CTR_input_condition,
            CTR_input_cluster,
            intracellular_network_condition,
            intracellular_network_cluster)

        return tf