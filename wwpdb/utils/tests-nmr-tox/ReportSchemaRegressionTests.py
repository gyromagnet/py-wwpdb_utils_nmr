"""A rejected source must still have a schema-valid diagnostic representation."""
import json
import unittest
from pathlib import Path

from jsonschema import Draft7Validator, ValidationError


class ReportSchemaRegressionTests(unittest.TestCase):

    def test_unidentified_content_is_explicitly_unknown(self):
        path = Path(__file__).parent.parent / 'tests-nmr' / 'json-schema' / 'nmr-data-processing-report-schema-v5.json'
        schema = json.loads(path.read_text(encoding='utf-8'))
        Draft7Validator.check_schema(schema)
        validator = Draft7Validator(schema['definitions']['input_source']['properties']['content_subtype'])
        validator.validate(None)
        validator.validate({'dist_restraint': 1})
        # Unknown is distinct from a known set of zero observations.
        validator.validate({})
        for invalid in ([], 'unknown', 0, {'dist_restraint': -1}):
            with self.subTest(invalid=invalid), self.assertRaises(ValidationError):
                validator.validate(invalid)
