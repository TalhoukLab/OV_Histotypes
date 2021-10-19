all: train_eval merge_eval iv_summary

sequential: seq_train_eval seq_merge_eval seq_iv_summary

# Train models and return evaluations
train_eval:
	./1-train_eval.sh $(filter-out $@,$(MAKECMDGOALS))

# Merge evaluations across reps and algorithms
merge_eval:
	./2-merge_eval.sh $(filter-out $@,$(MAKECMDGOALS))

# Internal validation summary and combined across datasets
iv_summary:
	./3-iv_summary.sh $(filter-out $@,$(MAKECMDGOALS))

# Sequential classifier
seq_train_eval:
	./sequential/1-seq_train_eval.sh $(filter-out $@,$(MAKECMDGOALS))

seq_merge_eval:
	./sequential/2-seq_merge_eval.sh $(filter-out $@,$(MAKECMDGOALS))

seq_iv_summary:
	./sequential/3-seq_iv_summary.sh $(filter-out $@,$(MAKECMDGOALS))

# Clean target
%:
		@:
clean:
	./assets/clean.sh
