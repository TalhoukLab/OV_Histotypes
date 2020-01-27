all: train_eval merge_eval iv_summary

# Train models and return evaluations
train_eval:
	./1-train_eval.sh $(filter-out $@,$(MAKECMDGOALS))

# Merge evaluations across reps and algorithms
merge_eval:
	./2-merge_eval.sh $(filter-out $@,$(MAKECMDGOALS))

# Internal validation summary and combined across datasets
iv_summary:
	./3-iv_summary.sh $(filter-out $@,$(MAKECMDGOALS))

# Clean target
%:
		@:
clean:
	./assets/clean.sh
