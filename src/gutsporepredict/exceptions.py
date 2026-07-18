"""Custom exceptions used by GutSporePredict."""


class GutSporePredictError(Exception):
    """Base exception for expected GutSporePredict failures."""


class ConfigurationError(GutSporePredictError):
    """Raised when configuration is invalid or incomplete."""


class InputValidationError(GutSporePredictError):
    """Raised when an input file or value is invalid."""


class ExternalToolError(GutSporePredictError):
    """Raised when an external bioinformatics tool fails."""


class OutputValidationError(GutSporePredictError):
    """Raised when expected output files are absent or invalid."""
