import React, { PureComponent } from "react";
import { connect } from "react-redux";
import { Drawer, H3, H1, UL } from "@blueprintjs/core";

@connect((state) => {
  return { singleValueCategories: state.metadata.singleValueCategories };
})
class InfoDrawer extends PureComponent {
  constructor(props) {
    super(props);

    this.state = {
      // CHANGE BACK TO FALSE BEFORE MERGE
      isOpen: true,
    };
  }

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
      singleValueCategories,
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
          <div style={{ margin: 24 }}>
            <H1>{datasetTitle}</H1>
            {singleValueCategories.length > 0 ? (
              <>
                <H3>Dataset Metadata</H3>
                <UL>
                  {Array.from(singleValueCategories).map((pair) => {
                    return <div key={pair[0]}>{`${pair[0]}: ${pair[1]}`}</div>;
                  })}
                </UL>
              </>
            ) : null}
            {aboutURL ? (
              <>
                <H3>More Info</H3>
                <a href={aboutURL} target="_blank" rel="noopener noreferrer">
                  {aboutURL}
                </a>
              </>
            ) : null}
          </div>
        </Drawer>
      </div>
    );
  }
}
export default InfoDrawer;
