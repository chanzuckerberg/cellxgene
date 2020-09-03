import React, { PureComponent } from "react";
import { connect, shallowEqual } from "react-redux";
import { Drawer, H3, H1, UL, Classes } from "@blueprintjs/core";
import Async from "react-async";
import {
  selectableCategoryNames,
  createCategorySummaryFromDfCol,
} from "../../util/stateManager/controlsHelpers";

@connect((state) => {
  return {
    annoMatrix: state.annoMatrix,
    schema: state.annoMatrix.schema,
    datasetTitle: state.config?.displayNames?.dataset ?? "",
    aboutURL: state.config?.links?.["about-dataset"],
  };
})
class InfoDrawer extends PureComponent {
  static watchAsync(props, prevProps) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  constructor(props) {
    super(props);

    this.state = {
      isOpen: false,
    };
  }

  fetchAsyncProps = async (props) => {
    const { schema } = props.watchProps;
    const { annoMatrix } = this.props;

    const allCategoryNames = selectableCategoryNames(schema).sort();

    const nonUserAnnoCategories = allCategoryNames.map((catName) => {
      const isUserAnno = schema?.annotations?.obsByName[catName]?.writable;
      if (!isUserAnno) return annoMatrix.fetch("obs", catName);
      return nonUserAnnoCategories.push(null);
    });
    const singleValueCategories = (
      await Promise.all(nonUserAnnoCategories)
    ).reduce((acc, categoryData, i) => {
      const catName = allCategoryNames[i];

      const column = categoryData.icol(0);
      const colSchema = schema.annotations.obsByName[catName];

      const categorySummary = createCategorySummaryFromDfCol(column, colSchema);

      const { numCategoryValues } = categorySummary;
      //  Add to the array if the category has only one value
      if (numCategoryValues === 1) {
        acc.set(catName, categorySummary.allCategoryValues[0]);
      }
      return acc;
    }, new Map());

    return { singleValueCategories };
  };

  handleClick = () => {
    this.setState((state) => {
      return { isOpen: !state.isOpen };
    });
  };

  handleKeyPress = (e) => {
    if (e.key === "Enter") {
      this.handleCategoryClick();
    }
  };

  render() {
    const { isOpen } = this.state;
    const {
      position,
      aboutURL,
      datasetTitle,
      title,
      children,
      schema,
    } = this.props;

    return (
      <div
        role="menuitem"
        tabIndex="0"
        onKeyPress={this.handleKeyPress}
        onClick={this.handleClick}
      >
        {children}
        <Drawer {...{ isOpen, position, title }}>
          <Async
            watchFn={InfoDrawer.watchAsync}
            promiseFn={this.fetchAsyncProps}
            watchProps={{ schema }}
          >
            <Async.Pending>
              <InfoFormat skeleton {...{ datasetTitle, aboutURL }} />
            </Async.Pending>
            <Async.Rejected>
              {(error) => {
                console.error(error);
                return <span>Failed to load info</span>;
              }}
            </Async.Rejected>
            <Async.Fulfilled>
              {(asyncProps) => {
                const { singleValueCategories } = asyncProps;
                return (
                  <InfoFormat
                    {...{ datasetTitle, aboutURL, singleValueCategories }}
                  />
                );
              }}
            </Async.Fulfilled>
          </Async>
        </Drawer>
      </div>
    );
  }
}

const InfoFormat = ({
  datasetTitle,
  singleValueCategories = new Map([
    ["cellxgene", "cellxgene"],
    ["cellxgenecellxgene", "cellxgene cellxgene"],
    ["cellxgene cellxgene", "cellxgene"],
    ["cell", "cellxgene cellxgene"],
    ["xgene", "cellxgene"],
    ["cellcellxgenecellxgene", "cellxgene cellxgene cell"],
    ["cellxgene cellxgene cell", "cellxgene"],
    ["cellx", "cellxgene cellxgene"],
  ]),
  aboutURL = "thisisabouthtelengthofaurl",
  skeleton = false,
}) => {
  return (
    <div style={{ margin: 24 }}>
      <H1 className={skeleton ? Classes.SKELETON : null}>{datasetTitle}</H1>
      {singleValueCategories.size > 0 ? (
        <>
          <H3 className={skeleton ? Classes.SKELETON : null}>
            Dataset Metadata
          </H3>
          <UL>
            {Array.from(singleValueCategories).map((pair) => {
              return (
                <li
                  className={skeleton ? Classes.SKELETON : null}
                  key={pair[0]}
                >{`${pair[0]}: ${pair[1]}`}</li>
              );
            })}
          </UL>
        </>
      ) : null}
      {aboutURL ? (
        <>
          <H3 className={skeleton ? Classes.SKELETON : null}>More Info</H3>
          <a
            className={skeleton ? Classes.SKELETON : null}
            href={aboutURL}
            target="_blank"
            rel="noopener noreferrer"
          >
            {aboutURL}
          </a>
        </>
      ) : null}
    </div>
  );
};
export default InfoDrawer;
