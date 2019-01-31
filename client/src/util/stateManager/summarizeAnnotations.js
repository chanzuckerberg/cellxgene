import _ from "lodash";
import finiteExtent from "../finiteExtent";

/*
Build and return obs/var summary using any annotation in the schema

Summary information for each annotation, keyed by annotation name.
Value will be an object, containing summary information.

For continuous annotations (int, float, etc):
  <annotation_name>: {
    categorical: false,
    range {
      min: <number>,
      max: <number>
    }
  }

For categorical annotations (boolean, string, category):
  <annotation_name>: {
    categorical: true,
    categories: [ <category1>, <category2>, ... ]
    categoryCounts: Map {
      <category1>: <number>,
      ...
    },
    numCategories: <number>
  }

Summarize will be returned for BOTH obs and var annotations.

Example:
  {
    "Splice_sites_Annotated": {
      categorical: false,
      range: {
        "min": 26,
        "max": 1075869
      }
    },
    "Selection": {
      categorical: true,
      numCategories, 3,
      categories: [ "Astrocytes(HEPACAM)", "Endothelial(BSC)", "Unpanned" ],
      categoryCounts: Map {
        "Astrocytes(HEPACAM)": 714,
        "Endothelial(BSC)": 123,
        "Unpanned": 665
      }
    }
  }

NOTE: will not summarize the required 'name' annotation, as that is
specified as unique per element.
*/
function _summarizeAnnotations(_schema, df) {
  const summary = _(_schema) // lodash wrapping: https://lodash.com/docs/4.17.11#lodash
    .filter(v => v.name !== "name") // don't summarize name
    .filter(v => !!df.col(v.name)) // ensure data contains this field
    .keyBy("name")
    .mapValues(anno => {
      const { name, type } = anno;
      const continuous = type === "int32" || type === "float32";
      const numRows = df.length;
      const col = df.col(name) ? df.col(name).asArray() : null;

      if (continuous) {
        let min;
        let max;
        let nan = 0;
        let pinf = 0;
        let ninf = 0;
        if (col) {
          for (let r = 0; r < numRows; r += 1) {
            const val = Number(col[r]);
            if (Number.isFinite(val)) {
              if (min === undefined) {
                min = val;
                max = val;
              } else {
                min = val < min ? val : min;
                max = val > max ? val : max;
              }
            } else if (Number.isNaN(val)) {
              nan += 1;
            } else if (val > 0) {
              pinf += 1;
            } else {
              ninf += 1;
            }
          }
        }
        return {
          categorical: false,
          range: { min, max, nan, pinf, ninf }
        };
      }

      /* else categorical */
      const categoryCounts = new Map();
      if (col) {
        for (let r = 0; r < numRows; r += 1) {
          const val = col[r];
          let curCount = categoryCounts.get(val);
          if (curCount === undefined) curCount = 0;
          categoryCounts.set(val, curCount + 1);
        }
      }
      return {
        categorical: true,
        categories: [...categoryCounts.keys()],
        categoryCounts,
        numCategories: categoryCounts.size
      };
    })
    .value();
  return summary;
}

export default function summarizeAnnotations(
  schema,
  obsAnnotations,
  varAnnotations
) {
  return {
    obs: _summarizeAnnotations(schema.annotations.obs, obsAnnotations),
    var: _summarizeAnnotations(schema.annotations.var, varAnnotations)
  };
}
