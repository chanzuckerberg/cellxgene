import { schema } from "./schema";
import { Dataframe, KeyIndex } from "../../../../src/util/dataframe";
import { encodeMatrixFBS } from "../../../../src/util/stateManager/matrix";

const indexedSchema = {
  obsByName: Object.fromEntries(
    schema.schema.annotations.obs.columns.map((v) => [v.name, v]) ?? []
  ),
  varByName: Object.fromEntries(
    schema.schema.annotations.var.columns.map((v) => [v.name, v]) ?? []
  ),
  embByName: Object.fromEntries(
    schema.schema.layout.obs.map((v) => [v.name, v]) ?? []
  ),
};

function makeMockColumn(s, length) {
  const { type } = s;
  switch (type) {
    case "int32":
      return new Int32Array(length).fill(Math.floor(99 * Math.random()));

    case "string":
      return new Array(length).fill("test");

    case "float32":
      return new Float32Array(length).fill(99 * Math.random());

    case "boolean":
      return new Array(length).fill(false);

    case "categorical":
      return new Array(length).fill(s.categories[0]);

    default:
      throw new Error("unkonwn type");
  }
}

function getEncodedDataframe(colNames, length, colSchemas) {
  const colIndex = new KeyIndex(colNames);
  const columns = colSchemas.map((s) => makeMockColumn(s, length));
  const df = new Dataframe([length, colNames.length], columns, null, colIndex);
  const body = encodeMatrixFBS(df);
  return body;
}

export function dataframeResponse(colNames, columns) {
  const colIndex = new KeyIndex(colNames);
  const df = new Dataframe(
    [columns[0].length, colNames.length],
    columns,
    null,
    colIndex
  );
  const body = encodeMatrixFBS(df);
  const headers = new Headers({
    "Content-Type": "application/octet-stream",
  });
  return () => Promise.resolve({ body, init: { status: 200, headers } });
}

function annotationObsResponse(request) {
  const url = new URL(request.url);
  const params = Array.from(url.searchParams.entries());
  const names = params
    .filter(([k]) => k === "annotation-name")
    .map(([, v]) => v);
  if (!names.every((n) => indexedSchema.obsByName[n])) {
    return Promise.reject(new Error("bad obs annotation name in URL"));
  }
  const colSchemas = names.map((n) => indexedSchema.obsByName[n]);
  const body = getEncodedDataframe(
    names,
    schema.schema.dataframe.nObs,
    colSchemas
  );

  const headers = new Headers({
    "Content-Type": "application/octet-stream",
  });
  return Promise.resolve({
    body,
    init: { status: 200, headers },
  });
}

function annotationVarResponse(request) {
  const url = new URL(request.url);
  const params = Array.from(url.searchParams.entries());
  const names = params
    .filter(([k]) => k === "annotation-name")
    .map(([, v]) => v);
  if (!names.every((n) => indexedSchema.varByName[n])) {
    return Promise.reject(new Error("bad var annotation name in URL"));
  }
  const colSchemas = names.map((n) => indexedSchema.varByName[n]);
  const body = getEncodedDataframe(
    names,
    schema.schema.dataframe.nVar,
    colSchemas
  );

  const headers = new Headers({
    "Content-Type": "application/octet-stream",
  });
  return Promise.resolve({
    body,
    init: { status: 200, headers },
  });
}

function layoutObsResponse(request) {
  const url = new URL(request.url);
  const params = Array.from(url.searchParams.entries());
  const names = params.filter(([k]) => k === "layout-name").map(([, v]) => v);
  if (!names.every((n) => indexedSchema.embByName[n])) {
    return Promise.reject(new Error("bad layout name in URL"));
  }
  const dims = names.map((n) => indexedSchema.embByName[n].dims).flat();
  const colSchemas = names
    .map((n) => [indexedSchema.embByName[n], indexedSchema.embByName[n]])
    .flat();
  const body = getEncodedDataframe(
    dims,
    schema.schema.dataframe.nObs,
    colSchemas
  );

  const headers = new Headers({
    "Content-Type": "application/octet-stream",
  });
  return Promise.resolve({
    body,
    init: { status: 200, headers },
  });
}

function dataVarResponse(request) {
  const url = new URL(request.url);
  const params = Array.from(url.searchParams.entries());

  const colNames = params.map((v) => `${v[0]}/${v[1]}`);
  const colSchemas = colNames.map(() => schema.schema.dataframe);
  const body = getEncodedDataframe(
    colNames,
    schema.schema.dataframe.nObs,
    colSchemas
  );

  const headers = new Headers({
    "Content-Type": "application/octet-stream",
  });
  return Promise.resolve({
    body,
    init: { status: 200, headers },
  });
}

export function responder(request) {
  const url = new URL(request.url);
  const { pathname } = url;
  if (pathname.endsWith("/annotations/obs")) {
    return annotationObsResponse(request);
  }
  if (pathname.endsWith("/annotations/var")) {
    return annotationVarResponse(request);
  }
  if (pathname.endsWith("/layout/obs")) {
    return layoutObsResponse(request);
  }
  if (pathname.endsWith("/data/var")) {
    return dataVarResponse(request);
  }
  return Promise.reject(new Error("bad URL"));
}

export function withExpected(expectedURL, expectedParams) {
  /*
  Do some additional error checking
  */
  return (request) => {
    // if URL is bogus, reject the promise
    const url = new URL(request.url);
    if (!url.pathname.endsWith(expectedURL)) {
      return Promise.reject(new Error("Unexpected URL!"));
    }
    const params = Array.from(url.searchParams.entries()).sort(
      (a, b) => a[0] < b[0]
    );
    expectedParams = expectedParams.slice().sort((a, b) => a[0] < b[0]);

    if (
      params.length !== expectedParams.length ||
      !params.every(
        (p, i) => p[0] === expectedParams[i][0] && p[1] === expectedParams[i][1]
      )
    ) {
      return Promise.reject(new Error("unexpected name requested in URL"));
    }

    return responder(request);
  };
}

export function annotationsObs(names) {
  return withExpected(
    "/annotations/obs",
    names.map((name) => ["annotation-name", name])
  );
}
