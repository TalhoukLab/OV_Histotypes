all: tune_wflows merge_results summarize_results

sequential: seq_tune seq_merge seq_summarize

retrain: retrain_tune retrain_merge retrain_summarize

gene_opt: gene_tune gene_merge gene_summarize

gene_opt_seq: gene_seq_tune gene_seq_merge gene_seq_summarize

tune_wflows:
	./1-tune_wflows.sh $(filter-out $@,$(MAKECMDGOALS))

merge_results:
	./2-merge_results.sh $(filter-out $@,$(MAKECMDGOALS))

summarize_results:
	./3-summarize_results.sh $(filter-out $@,$(MAKECMDGOALS))

# Sequential classifier
seq_tune:
	./sequential/1-tune_wflows.sh $(filter-out $@,$(MAKECMDGOALS))

seq_merge:
	./sequential/2-merge_results.sh $(filter-out $@,$(MAKECMDGOALS))

seq_summarize:
	./sequential/3-summarize_results.sh $(filter-out $@,$(MAKECMDGOALS))

# Retrain for sequential and two-step classifiers
retrain_tune:
	./retrain/1-tune_wflows.sh $(filter-out $@,$(MAKECMDGOALS))

retrain_merge:
	./retrain/2-merge_results.sh $(filter-out $@,$(MAKECMDGOALS))

retrain_summarize:
	./retrain/3-summarize_results.sh $(filter-out $@,$(MAKECMDGOALS))

# Gene optimization
gene_tune:
	./gene_opt/1a-tune_wflows.sh $(filter-out $@,$(MAKECMDGOALS))

gene_merge:
	./gene_opt/2a-merge_results.sh $(filter-out $@,$(MAKECMDGOALS))

gene_summarize:
	./gene_opt/3a-summarize_results.sh $(filter-out $@,$(MAKECMDGOALS))

# Gene optimization for sequential classifiers
gene_seq_tune:
	./gene_opt/1b-seq_tune_wflows.sh $(filter-out $@,$(MAKECMDGOALS))

gene_seq_merge:
	./gene_opt/2b-seq_merge_results.sh $(filter-out $@,$(MAKECMDGOALS))

gene_seq_summarize:
	./gene_opt/3b-seq_summarize_results.sh $(filter-out $@,$(MAKECMDGOALS))

# Clean target
%:
		@:
clean:
	./assets/clean.sh
