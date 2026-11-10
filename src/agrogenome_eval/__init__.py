
"""
AgroGenome-Eval

A scientific Python package for plant genomics,
genetic variant analysis, and AI evaluation.

Initial supported species:
- Zea mays (Maize)
- Triticum aestivum (Wheat)

Project status: Research and development
"""

from .variant_classifier import classify_variant

__all__ = ["classify_variant"]
