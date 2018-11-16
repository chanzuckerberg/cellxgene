# cellxgene REST API 0.2 specification

Items marked as (_future_) are intended for future implementation, and are included in the design to round out the concept, and highlight what we would do when/if we needed more functionality. The (_future_) items are not currently used by the cellxgene web application, and may be omitted from any backend - see [Current Front-End Dependencies](#current-front-end-dependencies) for more details.

_Caveat emptor, partial spec_: this is a sketch for a spec, not a full spec, and some shortcuts have been taken in the authorship. Best practices for a REST API are assumed but not documented here, such as API versioning, reasonable choices for HTTP response codes, etc. In addition, for clarity the JSON examples will not always have all required quoting (eg, on keys) - the actual implementation should use legal JSON/CSV.

_Note to readers_: please submit bugs, errata, and requests for clarifications as github issues on this repo.

## Terminology

Terminology mapping between the proposed model and other others:

- Data or dataframe - expression matrix
- Observation (**obs**) - cell, a column in Seurat
- Variable (**var**) - gene, a row in Seurat
- Annotation - metadata, which may be attached to either observations or variables, and will have a name
- Layout - N dimensional float coordinates (x,y,...) assigned to an obs or var
- Cluster - integer id assigned to a obs or var

Reverse map:

- Cell - observation or obs
- Gene - variable or var
- Metadata - annotation
- Expression - data
- Expression matrix - dataframe

## Underlying Data Model

The underlying data model is typical dataframe plus observation/variable annotations (aka metadata), and is more or less the same as the ScanPy [AnnData](http://anndata.readthedocs.io/en/latest/index.html) or Bioconductor dataframe model:

- Dataframe: expression matrix [nObs, nVar]
  - Observations are cells
  - Variables are genes (or equivalent variable, eg, isoforms)
  - Observations and variables are identified by an matrix coordinate pair (index, base 0)
  - Dataframe values must be defined (dense matrix)
  - The underlying type of the dataframe is IEEE single precision float (IEEE 754 "decimal32" data type, float32 in many language environments)
- Annotations (metadata):
  - Observations/cells or variables/genes can be annotated
  - Annotations are typed and described in a schema:
    - Obs or var annotation
    - Annotation name
    - Annotation type (all encoding is JSON format):
      - Float32 - scalar value, IEEE single precision float (decimal32), encoded as JSON number
      - Int32 - scalar value, 2's complement signed Int 32 bit encoded as JSON number
      - String - UTF-8 string, encoded as JSON string
      - Boolean - boolean value, encoded as JSON true/false
      - Categorical: value is an item from an _unordered_ enumeration. The enumeration schema is encoded as a JSON array. Enumeration categories may be a JSON number, string or boolean (or any combination thereof). Schema is encoded as a JSON array of enumeration category names.
  - Annotation values will conform to their schema, _or_ be undefined.
  - Annotation names are a string.
  - There is one guaranteed annotation name on every dimension: "name". This should be used as the default display name, must be a string and must be unique for the dimension (ie, all observation names are unique and all variable names are unique).

Missing from this model are ad hoc metadata/annotations not associated with a observation or variable - we can add that in the future if it is required (what ScanPy calls unstructured annotation).

### Assumptions about the Server

The following assumptions are made about server behavior:

- The server will have a data schema for annotations
- The server will have default algorithms for layout, clustering and differential expression (may add additional endpoints in the future)
- The server will serve a default data set, including annotations, cluster assignment and 2D layout (put another way, the client can request this data without any specification about how it is generated).

## Filters

There are several common filters supported by various routes. Route descriptions supporting these filters will reference the filter types.

### Annotation _Value_ Filters

Filter by the value of specific observation or variable annotation values. For example, the `/data` routes support sub-selection by the value of the annotations, which might be used to filter by some specific metadata (eg, `tissue_type: lung`).

Annotation value filters can be specified as a URL query parameter (using GET), or a JSON object in a request body (using PUT or POST)<sup>[1](#endnote-1)</sup>.

For a GET URL query parameter:

- Annotation name is encoded as `obs:name` or `var:name`<sup>[2](#endnote-2)</sup>.
- Enumerated values (string, categorical, boolean) are encoded as option lists, ie, `var:tissue=lung, obs:tumor=true`
- Scalar values (int32, float32) are encoded as ranges, ie, `obs:num_reads=1000,10000` where either min or max may be replaced with an asterisk to indicate a half-open range.
- Index filters are not allowed within GET URL query parameter filters
- Logically, filters are ANDed, except for repeated annotation names which are ORed. For example, `?X=A&X=B&Y=1` is evaluated as `((X==A or X==B) and Y==1)`

Example selection for _lung_ and _heart_ tissue with more than 1000 reads:

```
    GET /data?obs:tissue=lung&obs:tissue=heart&obs:num_reads=1000,*
```

For a PUT / POST request body, the filter encoding is a simple JSON object:

- Annotation names are embedded within either a obs or var object
- Enumerated values (eg, string, boolean categorical) are encoded as an array of options
- Scalar values are encoded as an object with a min and max field

Example selection for _lung_ and _heart_ tissue with more than 1000 reads:

```
    PUT or POST
    --
    {
      "filter": {
        "obs": {
          "annotation_value": [
            { "name": "tissue", "values": ["lung", "heart"] },
            { "name": "num_reads", "min": 1000 },
            { "name": "tumor", "values": true }
          ]
        }
      }
    }
```

### Obs or Var _Index_ Filters

Several PUT / POST routes allow filtering by observation or variable index, eg, access to a subset of the dataframe. Obs or var indexes may be specified individually or as a range.

For example, the following filter requests observations 1, 99, 1000-2000 from the dataframe (all variables are implicitly requested because no var filter is provided):

```
PUT or POST
--
{
  filter: {
    obs: {
      index: [
        1, 99, [1000, 2000]
      ]
    }
  }
}
```

### Complex Filters

Several routes support multiple filter types, which may be combined into a single PUT / POST body. Complex filters are additive.

Example:

```
PUT or POST
--
{
  "filter": {
    "obs": {
      "annotation_value": [
        { "name": "sex", "values": [ "F" ] }
      ],
      "index": [
        [0, 10000]
      ]
    },
    "var": {
      // ...
    }
  }
}
```

## Routes

Summary (details follow):

- `/config` - return application configuration
- `/schema` - data schema for the dataframe and dataframe annotations
- `/annotations/obs`, `/annotations/var` - fetch observation or variable annotation values, with ability to filter by annotation name (ie, fetch a subset of annotation data)
- `/data/obs`, `/data/var` - fetch dataframe values, with ability to subselect by annotation value and/or obs/var indices.
- `/cluster/obs`, `/cluster/var` - fetch cluster assignments, with ability to request re-clustering of a data subset
- `/layout/obs`, `/layout/var` - fetch layout information (N dimensional coordinates), with ability to request re-layout of a data subset
- `/diffexp` - fetch differential expression stats for a subset of data.
- `/saveLocal` - request that the server persist a copy of the caller-specified data to local storage, with caller-specified name.

All routes will be implemented with a prefix: <code>/api/v<strong><em>major.minor</em></strong></code>, where <em>major.minor</em> is a semver version number. Currently <em>version </em> is <em>0.2</em>, so, for example, the config route will be <code>/api/v0.2/config</code>. Prefixes are omitted in this document for clarity and brevity (eg, <code>/config</code> should be read as <code>/api/v0.2/config</code>).

The following routes are mandatory:

- `/config` - GET
- `/schema` - GET
- `/annotations/obs`, `/annotations/var` - GET and PUT
- `/data/obs`, `/data/var` - GET and PUT

All other routes are optional. Routes that are optional will advertise their presence or non-presence in the `/config` response.

### Common Error Codes

Routes that request potentially long-running (non-interactive speed) computation should provide threshold hints in /config response, and must return a 403 error immediately when a request exceeds, in the server's view, a reasonable interactive compute time.

Erroneously specified filters (eg, non-existent annotation name or other malformed filter) should return a 400 HTTP error.

### GET /config

Configuration information to assist in front-end adaptation to underlying engine, available functionality, interactive time limits, etc.

Information will include:

- API endpoint availability (eg, differential expression not available) and related limitations, including:
  - Feature enabled/disabled:
    - POST /cluster (re-clustering) supported: true, false
    - PUT /layout (re-layout) supported: true, false
    - POST /diffexp (differential expression calculation) available: true, false
  - Interactive compute speed hints; _approximate_ limitation on request size before it becomes non-interactive:
    - Re-clustering: if supported, this is an integer number of obs/vars at which the algorithm is non-interactive. Omitted if POST /cluster not supported
    - Re-layout: integer number of obs/vars at which the re-layout is non-interactive. Omitted if PUT /layout not supported
    - Diff Exp: integer number of obs/vars at which the differential expression compute is non-interactive. Omitted if POST /diffexp not supported
- Human readable information:
  - Data set name / title - for display purposes.
  - Engine name and version - for display purposes.
- System-wide configuration parameters, eg, configuration the end user has specified via a CLI or server configuration system.

**Response body:**

- `features` - encoded as array of features. All optional features must be included. Each feature will be identified by an HTTP method and path **prefix**, and will have availability status. If not specified as unavailable, the feature must be implemented.
  - Limitations are optional, and specified as an additional parameter on a feature, encoded as a number value with key `interactiveLimit`. By convention, this is the maximum number of inputs at which point the request is likely to return an error indicating non-interactive compute request. For example, this might be the maximum number of observations specified in the request to `PUT /layout/obs`.
  - mandatory routes do not need to be specified
- `displayNames` - names will include `engine` and `dataset` display names, as a JSON string. These should be relatively short strings (eg, suitable for window titles or equivalent)
- `parameters` - system configuration parameters, specified as a key/value pair in a simple object. See [Currently Defined Parameters](#currently-defined-parameters) for additional information.

```
GET /config
--
200
{
  "config": {
    "features": [
      // all /cluster/* paths not implemented
      { "method": "POST", "path": "/cluster/", "available": false },
      {
        "method": "PUT",
        "path": "/layout/obs",
        "available": true,
        "interactiveLimit": 10000
      },
      { "method": "PUT", "path": "/layout/var", "available": false }
    ],
    "displayNames": {
      "engine": "ScanPy version 1.33",
      "dataset": "/home/joe/mouse/blorth.csv"
    },
    "parameters": {
      "max-category-items": 1000,
      "verbose": false
      // name: value
    }
  }
}
```

### GET /schema

Schema definition for dataframe and annotations.

**URL Parameters:** none

**Response body:** encoded as a JSON object:

- Dataframe:
  - Dimensionality and type of dataframe, including:
    - nObs - number of observations
    - nVar - number of variables
    - type - dataframe element type (currently required to be float32)
- Dataframe annotations, both obs & var
  - name - encoded as a string
  - type - one of `float32`, `string`, `int32`, `boolean` or `categorical`
  - categories - meaningful only for `categorical` type. Specifies the complete list of category names.

Example:

```
GET /schema
--
200
{
  "schema": {
    "dataframe": {
      "nObs": 383,
      "nVar": 19944,
      "type": "float32"
    },
    "annotations": {
      "obs": [
        { "name": "name", "type": "string" },
        { "name": "tissue_type", "type": "string" },
        { "name": "num_reads", "type": "int32" },
        { "name": "sample_name", "type": "string" },
        {
          "name": "clusters",
          "type": "categorical",
          "categories"=[ 99, 1, "unknown cluster" ]
        },
        { "name": "QScore", "type": "float32" }
      ],
      "var": [
        { "name": "name", "type": "string" },
        { "name": "gene", "type": "string" }
      ]
    }
  }
}
```

### GET /annotations/obs, GET /annotations/var

Fetch annotations (metadata) for all observations or all variables. A subset of annotations may be selected by providing a query parameter `annotation-name` with a URL query list of annotation names.

Observations and variables are guaranteed to have a `name` annotation, which should contain a human-readable name for this element (eg. a gene or cell name). No other annotation names are required or may be assumed. The `name` annotation value must be unique for the dimension, ie, a given observation name is unique among all observation names and a given variable name is unique among all variables.

**URL Query Parameters:**

- annotations-name - URL query list of 1 or more annotation names. Only specified annotations are returned. If this query parameter is not specified, all annotations are returned.

**Response code:**

- 200 Success
- 400 Bad Request - one or more of the names specified with `annotation-name` are not associated with an annotation.

**Response body:** annotation description and values. Values conform to the schema returned by `/schema` and each record begins with the observation or variable index. Where the specific value is not defined, `null`will be returned. Values will be sorted by index.

Example:

```
GET /annotations/obs?annotation-name=tissue_type&annotation-name=sex&annotation-name=num_reads&annotation-name=clusters
--
200
{
  "names": [
    "tissue_type", "sex", "num_reads", "clusters"
  ],
  "data": [
    [ 0, "lung", "F", 39844, 99 ],
    [ 1, "heart", "M", 83, 1 ],
    [ 49, "spleen", null, 2, "unknown cluster" ],
    // [ obsOrVarIndex, value, value, value, value ],
    // ...
  ]
}
```

### PUT /annotations/obs, (_future_) PUT /annotations/var

Same as the `GET /annotations` routes, with additional ability to filter by observation or variable index and annotation value (complex filter)

**Fragment:** same as `GET /annotations/{obs|var}`

**URL Query Parameters:**

- annotations-name - URL query list of one or more annotation names. Only specified annotations are returned. If this query parameter is not specified, all annotations are returned.

**Request body:** Complex filters, including

- Annotation value filter
- Observation or (_future_) variable index filter

**Response code:**

- 200 Success
- 400 Bad Request - malformed filter or one or more of the annotation-name identifiers were not associated with an annotation name.

**Response body:**

- Same as `GET /annotations/{obs|var}`

### GET /data/obs, (future) GET /data/var

Get data from the dataframe (ie expression values), optionally sub-selected by annotation value filter. Supports JSON (default) and CSV responses, with annotation value filtering to enable basic deep linking. Response data must be organized according to path-specified primary dimension:

- observation-major (for `/data/obs`)
- variable-major (for `/data/var`)

**URL Query Parameters:**

- Annotation value filter
- Accept-type: specifies an acceptable media type for simple content type negotiation. Semantics are similar to HTTP Accept header, but included in query params to make it easy to cut & paste URL. If parameter is present, MUST take precedence over the value(s) specified in the Accept request header. If the service does not support the requested format, it must reply with a 406 - Not Acceptable error response. This parameter differs from an Accept header in that it only accepts one media type value and not a priority list (it is named accept-type and not accept because the behavior differs slightly). Acceptable values are the MIME types (same semantics as Accept-Type HTTP header):
  - application/json
  - text/csv

**Response code:**

- 200 - Success
- 406 - Not Acceptable (accept-type unacceptable MIME type)
- 400 - Bad Request - malformed filter

**Response body:** contents of the filtered dataframe, in primary-axis-major order, encoded as a dense matrix (eg, `/data/obs` will return data by observation). Data is preceded by a list of secondary axis indices, using the same index encoding as Index Filters (ie, an array of indices and/or index ranges). In JSON, this will be encoded as an array of indices with key obs or var; in CSV, this will be the header of the CSV response body. In both encodings, data will be sorted by index (both primary and secondary axis).

Example JSON response to a request which has filtered by annotation value:

```
GET /data/obs?accept-type=application/json&obs:tissue=lung
--
200
{
  "var": [ [ 0, 20000 ] ],
  "obs": [
    [ 1, 39483, 3902, 203, 0, 0, ..., 28 ],   // length == 20001
    // [ obsIndex, value1, value2, ..., value20000 ]
    // ...
  ]
}
```

### PUT /data/obs, (future) PUT /data/var

Same as `GET /data/{obs|var}`, but also accepts sub-selection by obs/var index filter. Accept may be specified as an HTTP header, and will support JSON and CSV. JSON will be the default served if no Accept header is supplied.

**URL Query Parameters:** none

**Request body:** Complex filters, including

- Annotation value filter
- Observation or (_future_) variable index filter

**Response code:**

- 200 Success
- 400 Bad Request - malformed filter
- 406 - Not Acceptable (accept-type missing or unacceptable MIME type)

**Response body:** same as `GET /data/{obs|var}`

### (future) POST /cluster/obs, (future) POST /cluster/var

Generate cluster assignments for the caller-specified subset of data, as indicated by the filter. This operation implicitly requests a re-clustering operation to be performed on the specified data.

If re-clustering is not supported by the server, must return an HTTP 501 response. If, in the view of the server, the request will exceed a reasonable interactive time period, must immediately return HTTP 403 error (error return _before_ attempting computation).

**URL Query Parameters:** none

**Request body:** complex filter:

- Annotation value filter
- Obs/var index filter

**Response code:**

- 200 - Success
- 400 Bad Request - malformed filter
- 403 Forbidden - non-interactive request
- 501 - Re-clustering is not implemented

**Response body:**

- For 200 Success, returns numeric cluster ID for each observation or variable, organized as an array of `[index, clusterId]`, an sorted by index:

  ```
  {
    "cluster": [
      [ 0, 24 ],
      // [ obsOrVarIndex, clusterID ],
      // ...
    ]
  }
  ```

* For 403 Forbidden, must return an explanation of refusal in the response body to allow differentiation from an authorization failure (see 403 response definition in [RFC 7231](https://tools.ietf.org/html/rfc7231#section-6.5.3))

### GET /layout/obs, (_future_) GET /layout/var

Get the _default_ layout for all observations or (_future_) all variables. Return dimensionality of layout and a coordinate list. The client must be prepared to accept dimensionality different than it natively displays, and should perform a best-effort attempt to display the higher dimensionality layout (eg, create a 2D view of a 3D layout).

**Response body:** dimensionality and coordinate list for each observation or variable, where coordinates are encoded as floats [0,1], sorted by index. Example 2D response:

```
{
  "layout": {
    "ndims": 2,
    "coordinates": [
      [ 0, 0.284483, 0.983744 ],
      [ 1, 0.038844, 0.739444 ],
      // [ obsOrVarIndex, X_coord, Y_coord ],
      // ...
    ]
  }
}
```

**Response code:**

- 200 - success
- 500 Internal Server Error - unprepared data, layout has not been stored in the input data

### (_future_) PUT /layout/obs, (_future_) PUT /layout/var

Generate layout for the caller-specified subset of data, as indicated by the filter. This operation implicitly requests a re-layout operation to be performed on the specified data. This operation will _commonly_ return the same results for any given caller-specified filter, but this behavior is not guaranteed.

If re-layout is not supported by the server, must return an HTTP 501 response. If, in the view of the server, the request will exceed a reasonable interactive time period, must immediately return HTTP 403 error (error return _before_ attempting computation).

**URL Query Parameters:** none

**Request Body:** complex filter

- Annotation value filter
- Observation or variable index filter

**Response code:**

- 200 - success
- 400 Bad Request - malformed filter
- 403 Forbidden - non-interactive request
- 501 - re-layout is not implemented

**Response body:**

- For 200 Success, same as `GET /layout/{obs|var}`
- For 403 Forbidden, must return an explanation of refusal in the response body to allow differentiation from an authorization failure (see 403 response definition in [RFC 7231](https://tools.ietf.org/html/rfc7231#section-6.5.3))

### POST /diffexp/obs

Generate differential expression (DE) statistics for two specified subsets of data, as indicated by the two provided observation complex filters (each filter provided either as annotation and/or obs/var index). **If only one subset of data is specified, the default assumption is that the second subset is all other observations (cells).** This operation implicitly requests a differential expression operation to be performed on the specified data.

Two modes are provided:

- `topN`: return top N differentially expressed variables (across all variables)
- `varFilter`: return DE for caller-provided variable filter (_future_)

Both modes perform calculations using a subset of observations, where each subset is defined by an observation filter (`set1` and `set2`). These filters must not include a variable filter.

If differential expression is not supported by the server, must return an HTTP 501 response. If, in the view of the server, the request will exceed a reasonable interactive time period, must immediately return HTTP 403 error (error return _before_ attempting computation).

**Request body:**

- Mode: `topN` or `varFilter` (specified with complex variable filter). `varFilter` must not include an "obs" section for the filter. Can be specified as one of:

  ```
  { "mode": "topN", "count": N }
  ```

  or

  ```
      {
        "mode": "varFilter",
        "varFilter": {
          "filter": {
            "var": {
              "annotation_value": [
                { "name": "gene", "values": [ "FOXP2", "ACTN1" ] }
              ]
            }
          }
        }
      }
  ```

- count: required for `topN` mode, and specifies the number of results to return.
- varFilter: required for `varFilter` mode, and specifies a variable filter.
- set1 - Complex filter - defines observation "set1".

  ```
  {
     "set1" {
         "filter": {
             "obs": {
                "annotation_value": [
                     { "name": "tissue", "values": [lung, heart] },
                     { "name": "num_reads", "min": 1000 },
                     { "name": "tumor", "values": true }
                ]
             }
         }
     }
  }
  ```

* set2 - Complex filter - defines observation "set2" (set 2 is optional, and defaults to all **other** observations, ie, those disjoint from set1)

**Response code:**

- 200 - success
- 400 Bad Request - malformed filter
- 403 Forbidden - non-interactive request
- 501 - differential expression is not implemented

**Response body:**

- For 200 Success, differential expression statistics returned as array of arrays, where each contains the following values:

  - **varIndex**: variable index for the computed results
  - **logfoldchange**: log fold-change of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group,
  - **pVal**: unadjusted p-value,
  - **pValAdj**: adjusted p-value

  Values ordered as:

  _varIndex_, _logfoldchange_, _pVal_, _pValAdj_

  For example:

  ```
  [
    [ 1720, 2.4679039, 2.3124478092035228e-175, 4.250279073316075e-172 ]
    // ...
  ]
  ```

- For 403 Forbidden, must return an explanation of refusal in the response body to allow differentiation from an authorization failure (see 403 response definition in [RFC 7231](https://tools.ietf.org/html/rfc7231#section-6.5.3))

Example:

```
POST /diffexp/obs
{
  "diffexp": {
    "set1": {
      "filter": {
        "obs": { "index": [ [0, 10000] ] }
      }
    },
    "set2": {
      "filter": {
        "obs": { "index": [28448, 4, [88888, 99999] ] }
      }
    },
    "mode": "topN",
    "count": 5
  }
}
---
200 - Success
{
  "diffexp": [
    [ 328, -2.569489, 2.655706e-63, 3.642036e-57 ],
    // [ varIdx, logfoldchange, pVal, pValAdj ],
    // ...
  ]
}
```

### (future) POST /data/saveSelection

Instruct the server to save or otherwise process the selected dataframe subset, as described the the filter. "Saved" data may optionally be given a name. If no name is provided, the server should pick a reasonable default. Action is synchronous - ie, a 200 return indicates completion of the request. Response will contain a human readable text message to inform the user the results of their action.

The intent of this route is to provide the user with a means of "saving" or otherwise processing filtered/selected data, for other uses, in a server-local resource. The precise semantics and mechanism are server-dependent, and the client should make no assumptions about server-side UI, data formats or persistence mechanism. Additionally, the client should make no assumptions about the lifespan or long-term persistence of the "saved" data.

For example, a server may provide:

- Store data subset in a CSV for sharing with other toolchains
- Persist a selection view within an HDF5 dataset region reference
- Pipe the selected data to another toolchain for immediate processing (ie, Unix pipe-like behavior)

(these are just examples - the server may implement any data export or integration UI, completely independent of client UI or client knowledge).

**URL Query Parameters:** name parameter, which for portability should only contain characters legal in a Posix file name (`A–Z a–z 0–9 . _ -`).

**Request body:** Complex filters, including

- Annotation value filter
- Observation or (_future_) variable index filter

**Response code:**

- 200 Success
- 400 Bad Request - malformed filter
- 500 - unable to perform requested action

**Response body:**

For 200 Success and 400 Bad Request:

```
{
  message: "Your results have been saved to /tmp/foo.csv"
}
```

Example:

```
POST /data/saveSelection?name=myFavCluster
{
  "filter": {
    "obs": {
      "annotation_value": [
        { "name": "cluster", "values": [ "99" ] }
      ],
    }
  }
}
--
200 Success
{
  "message": "OK"
}
```

## Current Front-End Dependencies

The cellxgene web front-end currently uses the following routes & features.

Routes:

- `GET /config`
- `GET /schema`
- `GET /annotations/obs`
- `GET /annotations/var`
- `GET /layout/obs` - get the default layout
- `PUT /data/obs` - request will contain a filter by var `name`
- `POST /diffexp/obs` - mode `topN`, typically with a `count` of 10, and two sets defined by an obs index filter (`{ filter: { obs: { index: [...] } } }`)

Requests include the following content negotiation headers:

- `Accept: application/json` - CSV not currently used in any routes
- `Accept-Encoding: gzip, deflate, br`

Currently unused routes:

- `/data/var`
- `/cluster/*`
- `/layout/var`
- `/data/saveSelection`

## Currently Defined Parameters

The system currently defines the following parameters:

- `max-category-items` - if any categorical annotation has a unique value count exceeding
  this parameter, it will not be displayed to the user. Number.
- `layout` - algorithm to use for graph layout. String.
- `diffexp` - algorithm to used to calculate differential expression. String.
- `title` - data set name. String.
- `debug` - run in debug mode. Boolean.
- `verbose` - more verbose logging. Boolean.

<!-- Endnotes themselves at the bottom. -->

## Notes

<a name="endnote-1">[1]</a>: Filters can often be large and complex, and exceed the URL length limitations - the POST method provides a means for the client to accommodate these cases. The GET method is provided for those cases where this is not concern, and/or where the client requires deep linking (subject to this URL length limit).

<a name="endnote-2">[2]</a>: Note to implementers -- this scheme requires careful use of URL encoding to ensure that both annotation names and annotation values are properly escaped. It requires that both normal URL encoding is applied, as well as ensuring that the annotation name has any embedded colon characters escaped.
