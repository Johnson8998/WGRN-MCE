%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data
clear
clc

load('TCGA_LUSC_Sample_data.mat')
load('TCGA_gene.mat')
load('LCC_human_gene.mat')
load('TCGA_LUSC_Methylation.mat') % Intersected with the samples of the expression data, and the order has been adjusted to be consistent.

% Read MCE1 transition probability and gene data
load('MCE1_Transition_Probability_TCGA_LUSC_normal.mat');
load('MCE1_gene_TCGA_LUSC_normal.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract normal samples based on SampleType
normal_sample_pos = strcmp(TCGA_LUSC_Methylation.SampleType,'normal'); % 8
% no_sample_pos = cellfun(@isempty,TCGA_LUSC_Methylation.SampleType); % 99

normal_samples = TCGA_LUSC_Methylation.Sample(normal_sample_pos);

normal_sample_pos1 = find(normal_sample_pos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read gene expression data
expression_data = TCGA_LUSC_Sample_data.ExpressionInOwnOrder;

a = TCGA_LUSC_Methylation.GeneID; % 15709/15740, there are duplicates
b = LCC_human_gene.GeneID; % 19892
c = TCGA_gene.GeneID(TCGA_gene.selected_POS_inLCC > 0); % 19129, all are in LCC_human_gene


% Ensure gene symbols match
unique_genes = b; 

% Initialize the MFE matrix
num_genes = length(unique_genes);
num_samples = length(normal_samples);
MFE_values = zeros(num_genes, num_samples);

% Traverse each sample
for sample_idx = 1:num_samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    geneID_inMCE = MCE1_gene_TCGA_LUSC_normal.GeneID(MCE1_gene_TCGA_LUSC_normal.POSinLCC(:,sample_idx));
    geneID_inTCGA = TCGA_gene.GeneID(TCGA_gene.selected_POS_inLCC > 0);

    gene_pos = ismember(geneID_inTCGA,geneID_inMCE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get the gene expression values of the current sample
    pi_values = expression_data(:, normal_sample_pos1(sample_idx));
    pi_values(~gene_pos) = 0;
    pi_values = pi_values ./ diag(sum(pi_values,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Traverse each gene
    for gene_idx = 1:num_genes   
        if ~(MCE1_gene_TCGA_LUSC_normal.POSinLCC(gene_idx,sample_idx))
            continue;  % This gene is not in the methylated network, skip it
        end

        % Initialize the MFE value
        mfe = 0;

        % Find the MCE transition matrix of this sample
        Matrix_coordinate = MCE1_Transition_Probability_TCGA_LUSC_normal.Matrix_coordinate{sample_idx};
        % Find the position of this gene in the MCE transition matrix
        gene_idx_1 = sum(MCE1_gene_TCGA_LUSC_normal.POSinLCC(1:gene_idx,sample_idx));
        % Find all non - zero neighbors of this gene
        neighbor = find(Matrix_coordinate(:,1) == gene_idx_1);
        % Get the weights (out - degrees) from this gene to its first - order neighbor genes
        pij = MCE1_Transition_Probability_TCGA_LUSC_normal.Matrix_Weight{sample_idx}(neighbor);


        % Find the position of this gene in the expressed genes TCGA_gene
        gene_idx_2 = find(TCGA_gene.selected_POS_inLCC(TCGA_gene.selected_POS_inLCC > 0) == gene_idx);


        % Calculate the weighted sum of expression values
        for neighbor_id = 1:length(neighbor)
            weight = pij(neighbor_id);  % Transition weight

            mfe = mfe + pi_values(gene_idx_2) * weight * log(pi_values(gene_idx_2) * weight);
        end
        
        % Store the MFE value of this gene
        MFE_values(gene_idx, sample_idx) = -mfe;
    end
end

% % % % % % Output as a table, with row names as gene symbols and column names as sample IDs
% % % % % output_table = array2table(MFE_values, 'VariableNames', normal_samples, 'RowNames', unique_genes);
% % % % % 
% % % % % % Save the results as a CSV file
% % % % % writetable(output_table, 'MFE_results.csv', 'WriteRowNames', true);