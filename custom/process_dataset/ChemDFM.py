import json
import random

type = "reaction_prediction"
file_path = f"custom/dataset/ChemDFM/{type}.jsonl"
result_path = f"custom/dataset/ChemDual/chemdfm_{type}.json"
instrctions = {
    "reaction_prediction": [
        # "Given the reactants and reagents below, come up with a possible product.",
        # "Please suggest a potential product based on the given reactants and reagents.",
        # "What product could potentially form from the reaction of the given reactants and reagents?",
        # "Given the following reactants and reagents, please provide a possible product.",
        # "Based on the given reactants and reagents, what product could potentially be produced?",
        # "Given the reactants and reagents listed, what could be a probable product of their reaction?",
        "Please provide a feasible product that could be formed using the given reactants and reagents.",
        # "Based on the given reactants and reagents, suggest a possible product.",
        # "Given the reactants and reagents provided, what is a possible product that can be formed?",
        # "Using the provided reactants and reagents, can you propose a likely product?",
        # "Using the listed reactants and reagents, offer a plausible product.",
        # "With the provided reactants and reagents, propose a potential product.",
    ],
    "retro_synthesis": [
        "What reactants could lead to the production of the following product?",
        "With the provided product, recommend some probable reactants that were likely used in its production.",
        "Can you identify some reactants that might result in the given product?",
        "Provide a list of potential reactants that may have produced the given product.",
        "Given the following product, please provide possible reactants.",
        "What are the possible reactants that could have formed the following product?",
        "Please suggest possible reactants for the given product.",
        "With the given product, suggest some likely reactants that were used in its synthesis.",
        "Provided the product below, propose some possible reactants that could have been used in the reaction.",
        "Given the product provided, propose some possible reactants that could have been employed in its formation.",
        "Given these product, can you propose the corresponding reactants?",
        "Please suggest potential reactants used in the synthesis of the provided product.",
        "Which reactants could have been used to generate the given product?",
        "Please suggest potential reactants for the given product.",
        "Based on the given product, provide some plausible reactants that might have been utilized to prepare it.",
    ],
}
data = []
with open(file_path) as f:
    for i in f.readlines():
        res = json.loads(i)
        data.append(
            {
                "instruction": random.choice(instrctions.get(type)),
                "input": res["input"],
                "output": ".".join(res["answer"]) if isinstance(res["answer"], list) else res["answer"],
            }
        )
with open(result_path, "w") as f:
    json.dump(data, f, indent=2)
