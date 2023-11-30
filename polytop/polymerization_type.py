class PolymerizationType:
    """
    Represents a pair of junction names specifying the polymerization type.

    Attributes:
    - junction_a (str): First junction name.
    - junction_b (str): Second junction name.
    """
    def __init__(self, junction_a: str, junction_b: str):
        """
        Initializes a PolymerJunction object with two junction names.
        """
        self.junction_a = junction_a
        self.junction_b = junction_b
