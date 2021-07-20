export const baseDataURL = "https://a.fake.url/api/v0.2";

(window as any).CELLXGENE = {
  API: {
    prefix: baseDataURL,
    version: "v0.2/",
  },
};

export { schema } from "./schema";
export * from "./routes";
